# Importing libraries 

library(torch)
library(diversitree)
library(luz)
library(ggplot2)
library(MLmetrics)
library(svMisc)
source("summary-statistics.R")
source("neural-network-functions.R")
source("infer-general-functions.R")


device <- "cuda:1"

# Parameters of phylogenetic trees
n_trees <- 10000 # total number of trees (train + valid + test)
n_taxa  <- 300 # size of the trees
lambda_range  <- c(0.1, 1.) # range of lambda values 
epsilon_range <- c(0.0, 0.9) # range of epsilon values 
ss_check <- TRUE

# Generate the trees and save 
out   <- load_dataset(n_trees, n_taxa, lambda_range, epsilon_range,
                      ss_check = ss_check)
trees <- out$trees # contains the phylogenetic trees generated 
vec.true.lambda <- out$lambda # contains the corresponding speciation rates 
vec.true.mu     <- out$mu # contains the corresponding extinction rates 

# Create the corresponding summary statistics data.frame
out       <- generate_ltt_dataframe(trees, n_taxa, vec.true.lambda, vec.true.mu)
df.ltt    <- out$ltt   # ltt dataframe 
df.rates  <- out$rates # rate dataframe
ds.ltt    <- convert_ltt_dataframe_to_dataset(df.ltt, df.rates)

# Parameters of the NN's training
n_train    <- 100
n_valid    <- 30
n_test     <- 20
n_epochs   <- 15 
batch_size <- 8

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

# Parameters of the RNN
n_hidden  <- 50  # number of neurons in hidden layers 
n_layer   <- 2   # number of stacked RNN layers 
p_dropout <- .01 # dropout probability


# Build the RNN 

rnn.net <- nn_module(
  initialize = function(n_input, n_hidden, n_layer, p_dropout = .01,
                        batch_first = TRUE) {
    self$rnn <- nn_lstm(input_size = n_input, hidden_size = n_hidden, 
                        dropout = p_dropout, num_layers = n_layer,
                        batch_first = batch_first)
    self$out <- nn_linear(n_hidden, 2)
  },
  
  forward = function(x) {
    x <- self$rnn(x)[[1]]
    x <- x[, dim(x)[2], ]
    x %>% self$out() 
  }
)

rnn <- rnn.net(1, n_hidden, n_layer, p_dropout) # create the RNN
rnn$to(device = device) # move the RNN to the choosen GPU 


net <- rnn.net(1, 50, 4, 0.)
net$to(device = device)


# Prepare training 

opt <- optim_adam(params = rnn$parameters) # optimizer 

train_batch <- function(b){
  opt$zero_grad()
  output <- rnn(b$x$reshape(c(b$x$shape, 1L))$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$backward()
  opt$step()
  loss$item()
}

valid_batch <- function(b) {
  output <- rnn(b$x$reshape(c(b$x$shape, 1L))$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$item()
}


# Training loop 

for (epoch in 1:n_epochs) {
  
  # Training part 
  rnn$train()
  train_loss <- c()
  
  coro::loop(for (b in train_dl) {
    loss <- train_batch(b)
    train_loss <- c(train_loss, loss)
  })
  
  cat(sprintf("epoch %0.3d/%0.3d - train - loss: %3.5f \n",
              epoch, n_epochs, mean(train_loss)))
  
  # Evaluation part 
  rnn$eval()
  valid_loss <- c()
  
  coro::loop(for (b in test_dl) {
    loss <- valid_batch(b)
    valid_loss <- c(valid_loss, loss)
  })
  
  cat(sprintf("epoch %0.3d/%0.3d - valid - loss: %3.5f \n", epoch, n_epochs, mean(valid_loss)))
}


# Evaluation of the predictions of the RNN w/ test set 

rnn$eval()
vec.pred.lambda <- c()
vec.pred.mu     <- c()
vec.true.lambda <- c()
vec.true.mu     <- c()

# Compute predictions 
coro::loop(for (b in test_dl) {
  out <- rnn(b$x$reshape(c(b$x$shape, 1L))$to(device = device))
  pred <- as.numeric(out$to(device = "cpu")) # move the tensor to CPU 
  true <- as.numeric(b$y)
  vec.pred.lambda <- c(vec.pred.lambda, pred[1])
  vec.pred.mu     <- c(vec.pred.mu, pred[2])
  vec.true.lambda <- c(vec.true.lambda, true[1])
  vec.true.mu     <- c(vec.true.mu, true[2])
})

# Prepare plot 
name.list <- list("lambda", "mu")
pred.list <- list(vec.pred.lambda, vec.pred.mu)
true.list <- list(vec.true.lambda, vec.true.mu)

dir <- "figures/rnn-ltt/"
fname <- create_predictions_plot_fname(n_trees, n_taxa, lambda_range, epsilon_range,
                                       dir, "nn", n_layer = n_layer,
                                       n_hidden = n_hidden, n_train = n_train)

fname.rnn <- paste(path, fname, "-rnn-ltt", sep = "")
fname.mle <- paste(path, fname, "-mle", sep="")

# Plot DNN predictions 
plot_pred_vs_true(pred.list, true.list, name.list, "rnn", save = FALSE, fname = fname.rnn)

# Plot MLE predictions for comparison 
trees <- trees[test_indices]
plot_mle_predictions(trees, vec.true.lambda, vec.true.mu, 
                     lambda_range, mu_range, 
                     type = "diversitree",
                     save = FALSE, fname = fname.mle)




