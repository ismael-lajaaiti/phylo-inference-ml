# ------------------------------------------------------------
# INFERRING MACROEVOLUTIONARY RATES WITH CNN AND TREE ENCODING
# ------------------------------------------------------------


# Importing Libraries and Sources

library(torch)
source("infer-general-functions.R")
source("neural-network-functions.R")



# Defining global parameters

device <- "cuda:1" # GPU where to run computations 
nn_type <- "cnn-encode" # type of the model: Convolutional Neural Network w/
                        # graph encoding 

n_trees <- 10000 # total number of trees (train + valid + test)
n_taxa  <- c(100, 1000) # size of the trees
param.range <- list("lambda"  = c(0.1,1.),
                    "epsilon" = c(0.,0.9)) # parameters range 

# Generate the trees and save 
out        <- load_dataset_trees(n_trees, n_taxa, param.range, ss_check = TRUE)
trees      <- out$trees # extract phylogenetic trees generated
true.param <- out$param # extract true values of the parameters 

# Create the corresponding encoding of the trees
#tensor.encode <- generate_encoding(trees, n_taxa)
tensor.encode <- readRDS("trees-dataset/ntrees-10000-ntaxa-100-1000-lambda-0.1-1-epsilon-0-0.9-sscheck-TRUE-encode.rds")
#saveRDS(tensor.encode %>% as.matrix(), "trees-dataset/ntrees-10000-ntaxa-100-1000-lambda-0.1-1-epsilon-0-0.9-sscheck-TRUE-encode.rds")
ds.encode     <- convert_encode_to_dataset(tensor.encode, true.param)

# Parameters of the NN's training
n_train    <- 9000
n_valid    <- 500
n_test     <- 500
n_epochs   <- 100
batch_size <- 64
patience   <- 10

# Creation of the train, valid and test dataset
train_indices     <- sample(1:n_trees, n_train)
not_train_indices <- setdiff(1:n_trees, train_indices)
valid_indices     <- sample(not_train_indices, n_valid)
test_indices      <- setdiff(not_train_indices, valid_indices)

extract_elements <- function(list_of_vectors, indices_to_extract){
  l <- as.list(do.call(cbind, list_of_vectors)[indices_to_extract,] 
               %>% as.data.frame())
  return(l)
}

train_ds <- ds.encode(tensor.encode[, train_indices],
                   extract_elements(true.param, train_indices))
valid_ds <- ds.encode(tensor.encode[, valid_indices],
                   extract_elements(true.param, valid_indices))
test_ds  <- ds.encode(tensor.encode[, test_indices],
                      extract_elements(true.param, test_indices))

# Creation of the dataloader 
train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
test_dl  <- test_ds  %>% dataloader(batch_size=1,          shuffle=FALSE)


n_hidden <- 5
n_layer  <- 3
ker_size <- 3
n_input  <- 2*max(n_taxa)-1
n_out    <- length(param.range)
p_dropout <- 0.05

# Build the CNN

cnn.net <- nn_module(
  
  "corr-cnn",
  
  initialize = function(n_input, n_out, n_hidden, n_layer, ker_size) {
    self$conv1 <- nn_conv1d(in_channels = 1, out_channels = n_hidden, kernel_size = ker_size)
    self$conv2 <- nn_conv1d(in_channels = n_hidden, out_channels = 2*n_hidden, kernel_size = ker_size)
    self$conv3 <- nn_conv1d(in_channels = 2*n_hidden, out_channels = 2*2*n_hidden, kernel_size = ker_size)
    n_flatten  <- compute_dim_ouput_flatten_cnn(n_input, n_layer, ker_size)
    self$fc1   <- nn_linear(in_features = n_flatten * (2*2*n_hidden), out_features = 100)
    self$fc2   <- nn_linear(in_features = 100, out_features = n_out)
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
      nnf_dropout(p = p_dropout) %>%
      nnf_relu() %>%
      
      self$fc2()
  }
)

cnn <- cnn.net(n_input, n_out, n_hidden, n_layer, ker_size) # create CNN
cnn$to(device = device) # Move it to the choosen GPU

# Prepare training 

opt <- optim_adam(params = cnn$parameters) # optimizer 

train_batch <- function(b){
  opt$zero_grad()
  output <- cnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$backward()
  opt$step()
  loss$item()
}

valid_batch <- function(b) {
  output <- cnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$item()
}


# Initialize parameters for the training loop 
epoch     <- 1
trigger   <- 0 
last_loss <- 100


# Training loop 

while (epoch < n_epochs & trigger < patience) {
  
  # Training part 
  cnn$train()
  train_loss <- c()
  
  coro::loop(for (b in train_dl) {
    loss <- train_batch(b)
    train_loss <- c(train_loss, loss)
  })
  
  cat(sprintf("epoch %0.3d/%0.3d - train - loss: %3.5f \n",
              epoch, n_epochs, mean(train_loss)))
  
  # Evaluation part 
  cnn$eval()
  valid_loss <- c()
  
  coro::loop(for (b in test_dl) {
    loss <- valid_batch(b)
    valid_loss <- c(valid_loss, loss)
  })
  
  current_loss <- mean(valid_loss)
  if (current_loss > last_loss){trigger <- trigger + 1}
  else{
    trigger   <- 0
    last_loss <- current_loss
  }
  
  cat(sprintf("epoch %0.3d/%0.3d - valid - loss: %3.5f \n", epoch, n_epochs, current_loss))
  
  epoch <- epoch + 1 
}


# Saving model for reproducibility 

#cat("\nSaving model...")
#fname.model <- get_model_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range, 
#                                   n_layer, n_hidden, n_train, ker_size)
#torch_save(cnn, fname.model)
#cat(paste("\n", fname.model, " saved.", sep = ""))
#cat("\nSaving model... Done.")

# Evaluation of the predictions of the RNN w/ test set 

cnn$eval()
nn.pred <- vector(mode = "list", length = n_out)
names(nn.pred) <- names(true.param)

# Compute predictions 
coro::loop(for (b in test_dl) {
  out <- cnn(b$x$to(device = device))
  pred <- as.numeric(out$to(device = "cpu")) # move the tensor to CPU 
  #true <- as.numeric(b$y)
  for (i in 1:n_out){nn.pred[[i]] <- c(nn.pred[[i]], pred[i])}
})

# Prepare plot 
name.param <- c("lambda", "mu")

# Save neural network predictions 
#save_predictions(pred.list, true.list, nn_type, n_trees, n_taxa,
#                 lambda_range, epsilon_range, n_test, n_layer, 
#                 n_hidden, n_train, ker_size)                                              


# Plot Predictions 
true.param.test <- as.list(as.data.frame(do.call(cbind, true.param))[test_indices,])
fname.mle <- get_mle_preds_save_name(n_trees, n_taxa, param.range, ss_check)
mle.pred <- readRDS(fname.mle)
mle.pred.test <- as.list(as.data.frame(do.call(cbind, mle.pred))[test_indices,])
pred.param.test <- list("mle" = mle.pred.test)
pred.param.test[[nn_type]] <- nn.pred
param.range.ajusted <- list("lambda" = c(0.1,1.), "mu" =c(0.,.9))

plot_pred_vs_true_all(pred.param.test, true.param.test, name.param, param.range.ajusted)

plot_bars_mle_vs_nn(pred.list.all, true.list, nn_type, name.list, save = TRUE, n_trees, n_taxa, 
                    lambda_range, epsilon_range, n_test, n_layer, n_hidden, n_train, ker_size)