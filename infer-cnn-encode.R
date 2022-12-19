# CNN with CBLV encoding

#### Dependecies ####
source("R/libraries.R")
source("R/neural-network-functions.R")
source("R/utils.R")

#### Set-up ####
gpu_id <- "cuda:1"
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

#### Prepare neural network ####
epoch_max <- 50
patience <- 2
n_hidden <- 8
n_layer <- 4
ker_size <- 10
n_input <- 3 * max(range_tree_size)
n_out <- 2
p_dropout <- 0.0

create_cnn <- nn_module(
    initialize = function() {
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

#### Training loop ####
cnn_trained <- train_neural_network(create_cnn, gpu_id, train_dl, valid_dl)

#### Evaluation on the test set ####
cnn$eval()
cnn_predictions <- vector(mode = "list", length = n_out)
names(cnn_predictions) <- names(true_param)
coro::loop(for (b in test_dl) {
    out <- cnn(b$x$to(device = gpu_id))
    pred <- as.numeric(out$to(device = "cpu"))
    for (i in 1:n_out) {
        cnn_predictions[[i]] <- c(cnn_predictions[[i]], pred[i])
    }
})

saveRDS(nn_predictions, str_c(dir, "predictions/cnn-cblv-predictions.rds"))
