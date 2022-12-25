# Infer speciation rate with GNN from phylogeny graph.

# Dependencies
using MLUtils
using Flux
using RData
using Graphs
using GraphNeuralNetworks

# Utility functions
"Graph COO encoding"
function graph_COO(edge_data)
    s = round.(Int, edge_data[:, 1])
    t = round.(Int, edge_data[:, 2])
    (s, t)
end

function format_node_data(node_data)
    transpose(Matrix(node_data))
end

function format_param(raw_param)
    lambda_vec = raw_param["lambda"]
    q_vec = raw_param["q"]
    reshape(hcat(lambda_vec, q_vec), (2, 1))
end

# Read data & create GNN input
input_data = load("data/test-node-edge-df.rds")
raw_param = load("data/test-data-for-jl-gnn.rds")["param"]
all_param = format_param(raw_param)
all_graph = GNNGraph[]
for (i, data) in enumerate(input_data)
    edge_data = data["edge"]
    node_data = data["node"]
    s, t = graph_COO(edge_data)
    ndata = format_node_data(node_data)
    param = all_param[i, :]
    g = GNNGraph(s, t, ndata=(; x=ndata), gdata=(; y=param))
    push!(all_graph, g)
end

# Neural network set-up & training
n_node_feature = size(all_graph[1].ndata.x, 1)
device = Flux.cpu # change to GPU if available
model = GNNChain(GCNConv(n_node_feature => 64),
    BatchNorm(64),
    x -> relu.(x),
    GCNConv(64 => 64, relu),
    GlobalPool(mean),
    Dense(64, 2)) |> device
ps = Flux.params(model)
opt = Adam(1.0f-4)
train_graphs, test_graphs = MLUtils.splitobs(all_graph, at=0.8)
train_loader = DataLoader(train_graphs, batchsize=10, shuffle=true, collate=true)
test_loader = DataLoader(test_graphs, batchsize=10, shuffle=false, collate=true)
loss(g::GNNGraph) = mean((vec(model(g, g.x)) - g.y) .^ 2)
loss(loader) = mean(loss(g |> device) for g in loader)
for epoch in 1:10
    for g in train_loader
        g = g |> device
        grad = gradient(() -> loss(g), ps)
        Flux.Optimise.update!(opt, ps, grad)
    end
    @info (; epoch, train_loss=loss(train_loader), test_loss=loss(test_loader))
end
