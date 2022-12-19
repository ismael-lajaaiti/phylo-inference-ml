#' Train the neural network
train_neural_network <- function(create_nn,
                                 gpu_id,
                                 train_dl,
                                 valid_dl,
                                 epoch_max = 50,
                                 patience = 3,
                                 verbose = TRUE) {
    # Set-up.
    nn <- create_nn()
    nn$to(device = gpu_id)
    opt <- torch::optim_adam(params = nn$parameters)
    epoch <- 1
    trigger <- 0
    last_loss <- Inf
    best_nn <- deepcopy_nn(create_nn, nn)
    while (epoch <= epoch_max && trigger < patience) {
        # Training.
        nn$train()
        train_loss <- c()
        coro::loop(for (b in train_dl) {
            loss <- train_step(nn, opt, b, gpu_id)
            train_loss <- c(train_loss, loss)
        })
        if (verbose) {
            cat(stringr::str_c(
                sprintf(
                    "Epoch %0.2d/%0.2d - Train loss: ", epoch, epoch_max
                ),
                format(mean(train_loss), scientific = TRUE, digits = 3),
                "\n"
            ))
        }
        # Validation.
        nn$eval()
        valid_loss <- c()
        coro::loop(for (b in valid_dl) {
            loss <- valid_step(nn, b, gpu_id)
            valid_loss <- c(valid_loss, loss)
        })
        current_loss <- mean(valid_loss)
        if (current_loss > last_loss) {
            trigger <- trigger + 1
        } else {
            trigger <- 0
            last_loss <- current_loss
            best_nn <- deepcopy_nn(create_nn, nn)
        }
        if (verbose) {
            cat(stringr::str_c(
                sprintf(
                    "Epoch %0.2d/%0.2d - Valid loss: ", epoch, epoch_max
                ),
                format(current_loss, scientific = TRUE, digits = 3),
                "\n"
            ))
        }
        epoch <- epoch + 1
    }
    best_nn$eval()
    best_nn
}

#' Realize the training step of the neural network
train_step <- function(nn, opt, batch, gpu_id) {
    opt$zero_grad()
    output <- nn(batch$x$to(device = gpu_id))
    target <- batch$y$to(device = gpu_id)
    loss <- torch::nnf_mse_loss(output, target)
    loss$backward()
    opt$step()
    loss$item()
}

#' Realize the validation step of the neural network
valid_step <- function(nn, batch, gpu_id) {
    output <- nn(batch$x$to(device = gpu_id))
    target <- batch$y$to(device = gpu_id)
    loss <- torch::nnf_mse_loss(output, target)
    loss$item()
}

#' Create an immutable copy of a given neural network
deepcopy_nn <- function(create_nn, nn_to_copy) {
    nn_copied <- create_nn()
    state_dict <- nn_to_copy$state_dict()
    nn_copied$load_state_dict(state_dict)
    nn_copied
}

#' Number of neurons of the flattening layer in a CNN
compute_dim_ouput_flatten_cnn <- function(n_input, n_layer, kernel_size = 2) {
    k <- kernel_size - 1
    for (i in 1:n_layer) {
        n_input <- as.integer((n_input - k) / 2)
    }
    n_input
}

#' Convert a matrix of CBLV phylogeny encoding to a torch dataset
convert_encode_to_dataset <- function(cblv_tensor, true_param) {
    torch::dataset(
        initialize = function(cblv_tensor, true_param) {
            self$x <- torch::torch_tensor(cblv_tensor) # input
            self$y <- torch::torch_tensor(do.call(cbind, true_param)) # target
        },
        .getitem = function(i) {
            list(x = self$x[, i]$unsqueeze(1), y = self$y[i, ])
        },
        .length = function() {
            self$y$size()[[1]]
        }
    )
}

#' Convert a dataframe of summary statistics to a torch dataset
convert_sumstat_to_dataset <- function(df) {
    torch::dataset(
        initialize = function(df, param_name, stat_excluded = c()) {
            df <- na.omit(df)
            x <- df[!(colnames(df) %in% c(param_name, stat_excluded))]
            x <- as.matrix(x)
            self$x <- torch::torch_tensor(x)
            y <- as.matrix(df[param_name])
            self$y <- torch::torch_tensor(y)
        },
        .getitem = function(i) {
            list(x = self$x[i, ], y = self$y[i, ])
        },
        .length = function() {
            self$y$size()[[1]]
        }
    )
}
