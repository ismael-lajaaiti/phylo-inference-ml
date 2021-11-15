# Importing Libraries and Sources

library(torch)
source("infer-general-functions.R")
source("neural-network-functions.R")


device <- "cuda:2" # GPU where to run computations 

nn_type <- "cnn-ltt" # type of the model: Convolutional Neural Network w/ LTT

# Parameters of phylogenetic trees
n_trees <- 10000 # total number of trees (train + valid + test)
n_taxa  <- c(100, 1000) # size of the trees
lambda_range  <- c(0.1, 1.) # range of lambda values 
epsilon_range <- c(0.0, 0.9) # range of epsilon values 
ss_check <- TRUE

# Generate the trees and save 
out   <- load_dataset_trees(n_trees, n_taxa, lambda_range, epsilon_range,
                            ss_check = ss_check)
trees           <- out$trees # contains the phylogenetic trees generated 
vec.true.lambda <- out$lambda # contains the corresponding speciation rates 
vec.true.mu     <- out$mu # contains the corresponding extinction rates 

# Create the corresponding summary statistics data.frame
out       <- generate_ltt_dataframe(trees, n_taxa, vec.true.lambda, vec.true.mu)
df.ltt    <- out$ltt   # ltt dataframe 
df.rates  <- out$rates # rate dataframe
ds.ltt    <- convert_ltt_dataframe_to_dataset_cnn(df.ltt, df.rates)

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
train_ds <- ds.ltt(df.ltt[, train_indices], df.rates[train_indices, ])
valid_ds <- ds.ltt(df.ltt[, valid_indices], df.rates[valid_indices, ])
test_ds  <- ds.ltt(df.ltt[, test_indices] , df.rates[test_indices, ])

# Creation of the dataloader 
train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
test_dl  <- test_ds  %>% dataloader(batch_size=1,          shuffle=FALSE)


n_hidden <- 10
n_layer  <- 3
ker_size <- 5

# Build the CNN

cnn.net <- nn_module(
  
  "corr-cnn",
  
  initialize = function(n_input, n_hidden, n_layer, ker_size) {
    self$conv1 <- nn_conv1d(in_channels = 1, out_channels = n_hidden, kernel_size = ker_size)
    self$conv2 <- nn_conv1d(in_channels = n_hidden, out_channels = 2*n_hidden, kernel_size = ker_size)
    self$conv3 <- nn_conv1d(in_channels = 2*n_hidden, out_channels = 2*2*n_hidden, kernel_size = ker_size)
    n_flatten <- compute_dim_ouput_flatten_cnn(n_input, n_layer, ker_size)
    self$fc1 <- nn_linear(in_features = n_flatten * (2*2*n_hidden), out_features = 100)
    self$fc2 <- nn_linear(in_features = 100, out_features = 2)
  },
  
  forward = function(x) {
    x %>% 
      self$conv1() %>%
      nnf_relu() %>%
      nnf_avg_pool1d(2) %>%

      self$conv2() %>%
      nnf_relu() %>%
      nnf_avg_pool1d(2) %>%

      self$conv3() %>%
      nnf_relu() %>%
      nnf_avg_pool1d(2) %>%

      torch_flatten(start_dim = 2) %>%
      self$fc1() %>%
      nnf_relu() %>%
      
      self$fc2()
  }
)

n_input <- max(n_taxa)
cnn <- cnn.net(n_input, n_hidden, n_layer, ker_size) # create CNN
cnn$to(device = device) # Move it to the choosen GPU

# Prepare training 

opt <- optim_adam(params = cnn$parameters) # optimizer 

train_batch <- function(b){
  opt$zero_grad()
  output <- cnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_l1_loss(output, target)
  loss$backward()
  opt$step()
  loss$item()
}

valid_batch <- function(b) {
  output <- cnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_l1_loss(output, target)
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
vec.pred.lambda <- c()
vec.pred.mu     <- c()
vec.true.lambda <- c()
vec.true.mu     <- c()

# Compute predictions 
coro::loop(for (b in test_dl) {
  out <- cnn(b$x$to(device = device))
  pred <- as.numeric(out$to(device = "cpu")) # move the tensor to CPU 
  true <- as.numeric(b$y)
  vec.pred.lambda <- c(vec.pred.lambda, pred[1])
  vec.pred.mu     <- c(vec.pred.mu, pred[2])
  vec.true.lambda <- c(vec.true.lambda, true[1])
  vec.true.mu     <- c(vec.true.mu, true[2])
})

# Prepare plot 
name.list <- list("lambda", "mu")
true.list <- list("lambda" = vec.true.lambda, "mu" = vec.true.mu)
pred.list <- list("lambda" = vec.pred.lambda, "mu" = vec.pred.mu)

# Save neural network predictions 
#save_predictions(pred.list, true.list, nn_type, n_trees, n_taxa,
#                 lambda_range, epsilon_range, n_test, n_layer, 
#                 n_hidden, n_train, ker_size)                                              


# Plot Predictions 
trees_test <- trees[test_indices] # test trees
plot_together_nn_mle_predictions(pred.list, true.list, trees_test, nn_type, n_trees, n_taxa, 
                                 lambda_range, epsilon_range, n_layer, n_hidden, n_train,
                                 ker_size = ker_size, save = TRUE)


pred.list.mle <- get_mle_preds(trees_test)
pred.list.all <- list()
pred.list.all[[nn_type]] <- pred.list
pred.list.all[["mle"]]   <- pred.list.mle

plot_bars_mle_vs_nn(pred.list.all, true.list, nn_type, name.list, save = TRUE, n_trees, n_taxa, 
                    lambda_range, epsilon_range, n_test, n_layer, n_hidden, n_train, ker_size)

