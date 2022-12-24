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
    best_nn$to(device = gpu_id)
    best_nn
}

#' Compute neural network parameter predictions on the test set.
get_predictions <- function(nn_trained,
                            test_dl,
                            true_param,
                            model_type,
                            gpu_id,
                            save = FALSE,
                            verbose = FALSE) {
    valid_model <- c("cnn-cblv", "dnn-ss", "cnn-ltt", "rnn-ltt", "gnn-phy")
    if (!(model_type %in% valid_model)) {
        error_msg <- stringr::str_c(valid_model, collapse = ", ")
        error_msg <- stringr::str_c("`model_type` should be in: ", error_msg)
        stop(error_msg)
    }
    n_out <- length(names(true_param))
    nn_predictions <- vector(mode = "list", length = n_out)
    names(nn_predictions) <- names(true_param)
    coro::loop(for (b in test_dl) {
        out <- nn_trained(b$x$to(device = gpu_id))
        pred <- as.numeric(out$to(device = "cpu"))
        for (i in 1:n_out) {
            nn_predictions[[i]] <- c(nn_predictions[[i]], pred[i])
        }
    })
    if (save) {
        fname <- stringr::str_c(
            "bisse-1e6-phylogenies/",
            "predictions/",
            model_type,
            "-predictions.rds"
        )
        if (file.exists(fname)) {
            error_msg <- stringr::str_c(
                "File ",
                fname,
                " already exists. Cannot overwrite.",
                sep = "'"
            )
            stop(error_msg)
        }
        saveRDS(nn_predictions, fname)
        if (verbose) {
            message(stringr::str_c("File", fname, "saved.", sep = " "))
        }
    }
    nn_predictions
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

#' Convert a dataframe of LTT to a torch dataset.
convert_ltt_to_dataset <- function(ltt_df, true_param, nn_type) {
    nn_type_low <- stringr::str_to_lower(nn_type)
    if (!(nn_type_low %in% c("cnn", "rnn"))) {
        stop('`nn_type` should be either "CNN" or "RNN".')
    }
    if (nn_type_low == "cnn") {
        ds <- convert_ltt_to_dataset_cnn(ltt_df, true_param)
    } else {
        ds <- convert_ltt_to_dataset_rnn(ltt_df, true_param)
    }
    ds
}

#' Internal sub-function to convert a dataframe of LTT
#' to a torch dataset for CNN.
convert_ltt_to_dataset_cnn <- function(ltt_df, true_param) {
    torch::dataset(
        initialize = function(ltt_df, true_param) {
            array_ltt <- torch::torch_tensor(as.matrix(ltt_df))
            self$x <- array_ltt
            self$y <- torch::torch_tensor(do.call(cbind, true_param))
        },
        .getitem = function(i) {
            list(x = self$x[, i]$unsqueeze(1), y = self$y[i, ])
        },
        .length = function() {
            self$y$size()[[1]]
        }
    )
}

#' Internal sub-function to convert a dataframe of LTT
#' to a torch dataset for RNN.
convert_ltt_to_dataset_rnn <- function(ltt_df, true_param) {
    torch::dataset(
        initialize = function(ltt_df, true_param) {
            tensor_ltt <- torch::torch_tensor(as.matrix(ltt_df))
            self$x <- tensor_ltt$unsqueeze(3)
            self$y <- torch::torch_tensor(do.call(cbind, true_param))
        },
        .getitem = function(i) {
            list(x = self$x[, i], y = self$y[i, ])
        },
        .length = function() {
            self$y$size()[[1]]
        }
    )
}
