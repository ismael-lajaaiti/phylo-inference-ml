# DNN with Summary Statistics

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
input_ss <- readRDS(str_c(dir, "formatted/sumstat/bisse-sumstat-all.rds"))
input_ss <- add_params_to_sumstat(input_ss, true_param)
ds <- convert_sumstat_to_dataset(input_ss)
indices_list <- readRDS(str_c(dir, "indices-1e6-phylogenies.rds"))
train_indices <- indices_list$train
valid_indices <- indices_list$valid
test_indices <- indices_list$test
train_ds <- ds(input_ss[train_indices, ], names(true_param))
valid_ds <- ds(input_ss[valid_indices, ], names(true_param))
test_ds <- ds(input_ss[test_indices, ], names(true_param))
batch_size <- 64
train_dl <- train_ds %>% dataloader(batch_size = batch_size, shuffle = TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size = batch_size, shuffle = FALSE)
test_dl <- test_ds %>% dataloader(batch_size = 1, shuffle = FALSE)

#### Prepare neural network ####
n_in <- length(train_ds[1]$x)
n_out <- 2
n_hidden <- 100
p_dropout <- 0.0
create_dnn <- nn_module(
    initialize = function() {
        self$fc1 <- nn_linear(in_features = n_in, out_features = n_hidden)
        self$fc2 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
        self$fc3 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
        self$fc4 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
        self$fc5 <- nn_linear(in_features = n_hidden, out_features = n_out)
    },
    forward = function(x) {
        x %>%
            self$fc1() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            self$fc2() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            self$fc3() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            self$fc4() %>%
            nnf_relu() %>%
            nnf_dropout(p = p_dropout) %>%
            self$fc5()
    }
)

#### Training loop ####
dnn_trained <- train_neural_network(create_dnn, gpu_id, train_dl, valid_dl)

#### Evaluation on the test set ####
dnn_trained$to(device = gpu_id)
nn_predictions <- vector(mode = "list", length = n_out)
names(nn_predictions) <- names(true_param)
coro::loop(for (b in test_dl) {
    out <- dnn_trained(b$x$to(device = gpu_id))
    pred <- as.numeric(out$to(device = "cpu"))
    for (i in 1:n_out) {
        nn_predictions[[i]] <- c(nn_predictions[[i]], pred[i])
    }
})

saveRDS(nn_predictions, str_c(dir, "predictions/dnn-ss-predictions.rds"))
