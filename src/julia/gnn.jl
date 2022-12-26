# Infer speciation rate with GNN from phylogeny graph.

# Dependencies
using Flux
using GraphNeuralNetworks
using Graphs
using MLUtils
using Printf
using RData

# Utility functions
"""
Graph COO encoding.
"""
function graph_COO(edge_data)
    s = round.(Int, edge_data[:, 1])
    t = round.(Int, edge_data[:, 2])
    (s, t)
end

"""
Convert the raw node dataframe into a matrix with the right dimensions.
"""
format_node_data(node_data) = transpose(Matrix(node_data))

"""
Convert the raw parameters into a matrix with the right dimensions.
"""
function format_param(raw_param)
    lambda_vec = raw_param["lambda"]
    q_vec = raw_param["q"]
    hcat(lambda_vec, q_vec)
end

"""
Create the train, valid and test datasets.
"""
function create_ds(all_data, all_param, idx)
    ds = GNNGraph[]
    for i in idx
        data = all_data[i]
        edge_data = data["edge"]
        node_data = data["node"]
        s, t = graph_COO(edge_data)
        ndata = format_node_data(node_data)
        param = all_param[i, :]
        g = GNNGraph(s, t; ndata = (; x = ndata), gdata = (; y = param))
        push!(ds, g)
    end
    MLUtils.batch(ds)
end

# Read data & create GNN input
input_data = load("data/test-node-edge-df.rds")
raw_param = load("data/test-data-for-jl-gnn.rds")["param"]
all_param = format_param(raw_param)
all_graph = GNNGraph[]
indices = load("data/some-indices.rds")
ds_dict = Dict{}()
for split in ["train", "test", "valid"]
    idx = indices[split] |> collect # collect might not be necessary, to check
    ds_dict[split] = create_ds(input_data, all_param, idx)
end
train_dl = DataLoader(ds_dict["train"]; batchsize = 40, shuffle = true)
valid_dl = DataLoader(ds_dict["valid"]; batchsize = 10, shuffle = true)
test_dl = DataLoader(ds_dict["test"]; batchsize = 1, shuffle = false)

# Neural network set-up
n_node_feature = size(input_data[1]["node"], 2)
device = Flux.cpu # change to GPU if available
model =
    GNNChain(
        GCNConv(n_node_feature => 64),
        x -> relu.(x),
        GCNConv(64 => 64, relu),
        GlobalPool(mean),
        Dense(64, 2),
    ) |> device
ps = Flux.params(model)
opt = Adam(1.0f-4)
loss(g::GNNGraph) = mean((model(g, g.x) - g.y) .^ 2) # MSE
loss(loader) = mean(loss(g |> device) for g in loader)

# Neural network training
epoch_max = 10
min_valid_loss = Inf
patience = 3 # delay for early stopping
patience_count = 0
for epoch in 1:epoch_max
    # Training step
    for g in train_dl
        g = g |> device
        grad = gradient(() -> loss(g), ps)
        Flux.Optimise.update!(opt, ps, grad)
    end
    train_loss = loss(train_dl)
    train_loss_str = @sprintf "%.2E" train_loss

    # Validation step
    valid_loss = loss(valid_dl)
    valid_loss_str = @sprintf "%.2E" valid_loss

    # Early stopping check
    if valid_loss < min_valid_loss
        min_valid_loss = valid_loss
        patience_count = 0
    else
        patience_count += 1
    end
    epoch_str = @sprintf "%0.2d" epoch
    @info "Epoch = $epoch - Loss = (train = $train_loss_str, valid = $valid_loss_str) \
        - Patience = $patience_count"
    (patience_count >= patience) && break
end
