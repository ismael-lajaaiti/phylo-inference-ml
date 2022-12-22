# CNN with LTT coordinates

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
input_ltt <- readRDS(str_c(dir, "formatted/ltt/bisse-ltt-all.rds"))
ds <- convert_ltt_to_dataset(input_ltt, true_param, "RNN")
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
batch_size <- 32
train_dl <- train_ds %>% dataloader(batch_size = batch_size, shuffle = TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size = batch_size, shuffle = FALSE)
test_dl <- test_ds %>% dataloader(batch_size = 1, shuffle = FALSE)

#### Prepare the NN ####
n_hidden <- 100
n_layer <- 1
p_dropout <- 0.0
n_out <- 2


create_rnn <- nn_module(
    initialize = function() {
        self$rnn <- nn_lstm(
            input_size = 1,
            hidden_size = n_hidden,
            dropout = p_dropout,
            num_layers = n_layer,
            batch_first = TRUE
        )
        self$out <- nn_linear(max(range_tree_size) * n_hidden, n_out)
    },
    forward = function(x) {
        # Utility function.
        extract_output <- function(x) x[[1]]

        # Forward flow.
        x %>%
            self$rnn() %>%
            extract_output() %>%
            torch_flatten(start_dim = 2) %>%
            self$out()
    }
)

#### Training loop ####
rnn_trained <- train_neural_network(create_rnn, gpu_id, train_dl, valid_dl)

#### Evaluation on the test set ####
rnn_predictions <- get_predictions(
    rnn_trained,
    test_dl,
    true_param,
    "rnn-ltt",
    gpu_id,
    save = TRUE,
    verbose = TRUE
)
