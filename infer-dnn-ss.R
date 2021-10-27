# Importing libraries 

library(torch)
library(luz)
library(ggplot2)
library(MLmetrics)
source("summary-statistics.R")
source("infer-general-functions.R")
source("neural-network-functions.R")


# Parameters of phylogenetic trees
n_trees <- 10500 # total number of trees (train + valid + test)
n_taxa  <- 300 # size of the trees
lambda_range <- c(0.15, 0.25) # range within random lambda will be generated 
mu_range     <- c(0.0 , 0.1) # same for mu

# Generate the trees and save 
out   <- generate_trees(n_trees, n_taxa, lambda_range, mu_range, ss_check = TRUE)
trees <- out$trees # contains the phylogenetic trees generated 
vec.true.lambda <- out$lambda # contains the corresponding speciation rates 
vec.true.mu     <- out$mu # contains the corresponding extinction rates 

# Create the corresponding summary statistics data.frame
df <- generate_ss_dataframe_from_trees(trees, vec.true.lambda, vec.true.mu)

# Parameters of the NN's training
n_train    <- 9000
n_valid    <- 1000
n_test     <- 500
batch_size <- 64

# Creation of the train, valid and test dataset
ds <- convert_ss_dataframe_to_dataset(df)
train_indices <- sample(1:nrow(df), n_train)
not_train_indices <- setdiff(1:nrow(df), train_indices)
valid_indices <- sample(not_train_indices, n_valid)
test_indices  <- setdiff(not_train_indices, valid_indices)
train_ds <- ds(df[train_indices, ])
valid_ds <- ds(df[valid_indices, ])
test_ds  <- ds(df[test_indices, ])

# Creation of the dataloader 
train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
test_dl  <- test_ds  %>% dataloader(batch_size=batch_size, shuffle=FALSE)

# DNN parameters 
n_in      <- ncol(df) - 2 # number of neurons of the input layer 
n_hidden  <- 50 # number of neurons in the hidden layers 
p_dropout <- 0.01 # dropout probability 
n_epochs  <- 100 # maximum number of epochs for the training 
patience  <- 7 # patience of the early stopping 

# Build the DNN 
dnn <- nn_module(
  
  "ss-dnn", 
  
  initialize = function(){
    self$fc1 <- nn_linear(in_features = n_in, out_features = n_hidden)
    self$fc2 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
    self$fc3 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
    self$fc4 <- nn_linear(in_features = n_hidden, out_features = 2)
  }, 
  
  forward = function(x){
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
      
      self$fc4()
  }
)

cat("\nTraining DNN...")

# Fit the DNN
dnn.fit <- dnn %>%
  setup(
    loss = function(y_hat, y_true) nnf_mse_loss(y_hat, y_true),
    optimizer = optim_adam
  ) %>%
  fit(train_dl, epochs = n_epochs, valid_data = valid_dl, 
      callbacks = list(luz_callback_early_stopping(patience = patience)))

cat("\nTraining DNN... Done.")

# Evaluate DNN performance w/ predictions 
preds           <- predict(dnn.fit, test_dl)$to(device = "cpu") # get DNN preds
vec.pred.lambda <- as.matrix(preds)[,1] # extract lambda predictions 
vec.pred.mu     <- as.matrix(preds)[,2] # extract mu predictions 
vec.true.lambda <- df$lambda[test_indices] # corresponding lambda true values 
vec.true.mu     <- df$mu[test_indices] # corresponding mu true values 

# Prepare plots 
name.list <- list("lambda", "mu")
true.list <- list(vec.true.lambda, vec.true.mu)
pred.list <- list(vec.pred.lambda, vec.pred.mu)

n_layers <- length(dnn()$parameters)/2 - 2 
path <- "figures/dnn-ss/"
fname <- paste("ntaxa", n_taxa,
               "lambda", lambda_range[1], lambda_range[2], 
               "mu", mu_range[1], mu_range[2],
               "ntest", n_test,
               "nlayer", n_layers,
               "nhidden", n_hidden,
               "ntrain", n_train, sep="-")

fname.dnn <- paste(path, fname, "-dnn-ss", sep = "")
fname.mle <- paste(path, fname, "-mle", sep="")

# Plot DNN predictions 
plot_pred_vs_true(pred.list, true.list, name.list, "dnn", save = TRUE, fname = fname.dnn)

# Plot MLE predictions for comparison 
trees <- trees[test_indices]
plot_mle_predictions(trees, vec.true.lambda, vec.true.mu, 
                     lambda_range, mu_range, 
                     type = "diversitree",
                     save = TRUE, fname = fname.mle)
