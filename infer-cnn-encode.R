# CNN with CBLV encoding

#### Dependecies ####
source("R/libraries.R")
source("R/neural-network-functions.R")
source("R/utils.R")

#### Set-up ####
gpu_id <- "cuda:1"
nn_type <- "cnn-encode"
diversification <- "bisse"
range_tree_size <- c(100, 1000)

#### Read and prepare input data ####
dir <- "bisse-1e6-phylogenies/"
true_param <- readRDS(str_c(dir, "raw/bisse-true-params-all.rds"))
input_cblv <- readRDS(str_c(dir, "formatted/cblv/bisse-cblv-all.rds"))
ds <- convert_encode_to_dataset(input_cblv, true_param)
indices_list <- readRDS(str_c(dir, "indices-1e6-phylogenies.rds"))
train_indices <- indices_list$train
valid_indices <- indices_list$valid
test_indices <- indices_list$test

train_ds <- ds(
    input_cblv[, train_indices],
    extract_elements(true_param, train_indices)
)
valid_ds <- ds(
    input_cblv[, valid_indices],
    extract_elements(true_param, valid_indices)
)
test_ds <- ds(
    input_cblv[, test_indices],
    extract_elements(true_param, test_indices)
)

batch_size <- 64
train_dl <- train_ds %>% dataloader(batch_size = batch_size, shuffle = TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size = batch_size, shuffle = FALSE)
test_dl <- test_ds %>% dataloader(batch_size = 1, shuffle = FALSE)

# Prepare neural network.
epoch_max <- 50
patience <- 2
n_hidden <- 8
n_layer <- 4
ker_size <- 10
n_input <- 3 * max(range_tree_size)
n_out <- 2
p_dropout <- 0.0

create_cnn <- nn_module(
    "corr-cnn",
    initialize = function(n_input, n_out, n_hidden, n_layer, ker_size) {
        self$conv1 <- nn_conv1d(
            in_channels = 1,
            out_channels = n_hidden,
            kernel_size = ker_size
        )
        self$conv2 <- nn_conv1d(
            in_channels = n_hidden,
            out_channels = 2 * n_hidden,
            kernel_size = ker_size
        )
        self$conv3 <- nn_conv1d(
            in_channels = 2 * n_hidden,
            out_channels = 4 * n_hidden,
            kernel_size = ker_size
        )
        self$conv4 <- nn_conv1d(
            in_channels = 4 * n_hidden,
            out_channels = 8 * n_hidden,
            kernel_size = ker_size
        )
        n_flatten <- compute_dim_ouput_flatten_cnn(n_input, n_layer, ker_size)
        self$fc1 <- nn_linear(
            in_features = n_flatten * (8 * n_hidden),
            out_features = 100
        )
        self$fc2 <- nn_linear(
            in_features = 100,
            out_features = n_out
        )
    },
    forward = function(x) {
        x %>%
            self$conv1() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            nnf_avg_pool1d(2) %>%
            self$conv2() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            nnf_avg_pool1d(2) %>%
            self$conv3() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            nnf_avg_pool1d(2) %>%
            self$conv4() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            nnf_avg_pool1d(2) %>%
            torch_flatten(start_dim = 2) %>%
            self$fc1() %>%
            nnf_dropout(p = p_dropout) %>%
            nnf_relu() %>%
            self$fc2()
    }
)

cnn <- create_cnn(n_input, n_out, n_hidden, n_layer, ker_size)
cnn$to(device = gpu_id)
opt <- optim_adam(params = cnn$parameters) # optimizer

#### Training loop ####
epoch <- 1
trigger <- 0
last_loss <- Inf
best_cnn <- cnn

while (epoch < epoch_max && trigger < patience) {
    # Training.
    cnn$train()
    train_loss <- c()
    coro::loop(for (b in train_dl) {
        loss <- train_step(cnn, opt, b, gpu_id)
        train_loss <- c(train_loss, loss)
    })
    cat(str_c(
        sprintf(
            "Epoch %0.2d/%0.2d - Train loss: ", epoch, epoch_max
        ),
        format(mean(train_loss), scientific = TRUE, digits = 3),
        "\n"
    ))

    # Validation.
    cnn$eval()
    valid_loss <- c()
    coro::loop(for (b in test_dl) {
        loss <- valid_step(cnn, b, gpu_id)
        valid_loss <- c(valid_loss, loss)
    })
    current_loss <- mean(valid_loss)
    if (current_loss > last_loss) {
        trigger <- trigger + 1
    } else {
        trigger <- 0
        last_loss <- current_loss
        best_cnn <- cnn
    }
    cat(str_c(
        sprintf(
            "Epoch %0.2d/%0.2d - Valid loss: ", epoch, epoch_max
        ),
        format(current_loss, scientific = TRUE, digits = 3),
        "\n"
    ))
    epoch <- epoch + 1
}

#### Evaluation on the test set ####
best_cnn$eval()
nn_predictions <- vector(mode = "list", length = n_out)
names(nn_predictions) <- names(true_param)
coro::loop(for (b in test_dl) {
    out <- best_cnn(b$x$to(device = gpu_id))
    pred <- as.numeric(out$to(device = "cpu"))
    for (i in 1:n_out) {
        nn_predictions[[i]] <- c(nn_predictions[[i]], pred[i])
    }
})

saveRDS(str_c(dir, "predictions/cnn-cblv-predictions.rds"))
