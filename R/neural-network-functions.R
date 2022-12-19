#' Realize the training step of the neural network
train_step <- function(nn, opt, batch, gpu_id) {
    opt$zero_grad()
    output <- nn(batch$x$to(device = gpu_id))
    target <- batch$y$to(device = gpu_id)
    loss <- nnf_mse_loss(output, target)
    loss$backward()
    opt$step()
    loss$item()
}

#' Realize the validation step of the neural network
valid_step <- function(nn, batch, gpu_id) {
    output <- nn(batch$x$to(device = gpu_id))
    target <- batch$y$to(device = gpu_id)
    loss <- nnf_mse_loss(output, target)
    loss$item()
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
