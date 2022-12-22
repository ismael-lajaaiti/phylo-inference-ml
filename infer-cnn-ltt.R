# CNN with LTT coordinates

#### Dependecies ####
source("R/libraries.R")
source("R/neural-network-functions.R")
source("R/utils.R")

#### Set-up ####
gpu_id <- "cuda:3"
diversification <- "bisse"
range_tree_size <- c(100, 1000)

#### Read and prepare input data ####
dir <- "bisse-1e6-phylogenies/"
true_param <- readRDS(str_c(dir, "raw/bisse-true-params-all.rds"))
input_ltt <- readRDS(str_c(dir, "formatted/ltt/bisse-ltt-all.rds"))
ds <- convert_ltt_to_dataset(input_ltt, true_param, "CNN")
indices_list <- readRDS(str_c(dir, "indices-1e6-phylogenies.rds"))
train_indices <- indices_list$train
valid_indices <- indices_list$valid
test_indices <- indices_list$test
sliced_true_param <- lapply(
    indices_list,
    function(x) extract_elements(true_param, x)
)
train_ds <- ds(input_ltt[, train_indices], sliced_true_param$train)
valid_ds <- ds(input_ltt[, valid_indices], sliced_true_param$valid)
test_ds <- ds(input_ltt[, test_indices], sliced_true_param$test)
batch_size <- 64
train_dl <- train_ds %>% dataloader(batch_size = batch_size, shuffle = TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size = batch_size, shuffle = FALSE)
test_dl <- test_ds %>% dataloader(batch_size = 1, shuffle = FALSE)

#### Prepare neural network ####
n_hidden <- 8
n_layer <- 3
ker_size <- 5
p_dropout <- 0.0
n_input <- max(range_tree_size)
n_out <- 2
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
            out_channels = 2 * 2 * n_hidden,
            kernel_size = ker_size
        )
        n_flatten <- compute_dim_ouput_flatten_cnn(
            n_input,
            n_layer,
            ker_size
        )
        self$fc1 <- nn_linear(
            in_features = n_flatten * (2 * 2 * n_hidden),
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
            torch_flatten(start_dim = 2) %>%
            self$fc1() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            self$fc2()
    }
)

#### Training loop ####
cnn_trained <- train_neural_network(create_cnn, gpu_id, train_dl, valid_dl)

#### Evaluation on the test set ####
cnn_predictions <- get_predictions(
    cnn_trained,
    test_dl,
    true_param,
    "cnn-ltt",
    gpu_id,
    save = TRUE,
    verbose = TRUE
)
