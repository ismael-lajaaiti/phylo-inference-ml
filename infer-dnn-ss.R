# Importing libraries 

library(torch)
library(luz)
source("summary-statistics.R")
source("infer-general-functions.R")
source("neural-network-functions.R")


nn_type <- "dnn-ss" # type of the model: Deep Neural Network w/ Summary Statistics

# Parameters of phylogenetic trees
n_trees <- 10000 # total number of trees (train + valid + test)
n_taxa  <- c(100,1000) # size of the trees
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
df <- load_dataset_summary_statistics(n_trees, n_taxa, lambda_range,
                                      epsilon_range)

# Parameters of the NN's training
n_train    <- 9000
n_valid    <- 500
n_test     <- 500
batch_size <- 32

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
patience  <- 10 # patience of the early stopping 

# Build the DNN 
dnn <- nn_module(
  
  "ss-dnn", 
  
  initialize = function(){
    self$fc1 <- nn_linear(in_features = n_in, out_features = n_hidden)
    self$fc2 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
    self$fc3 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
    self$fc4 <- nn_linear(in_features = n_hidden, out_features = n_hidden)
    self$fc5 <- nn_linear(in_features = n_hidden, out_features = 2)
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
      
      self$fc4() %>%
      nnf_relu() %>%
      nnf_dropout(p = p_dropout) %>%
      
      self$fc5()
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

cat("\nSaving model...")
n_layer <- length(dnn()$parameters)/2 - 2 
# epochs <- length(dnn.fit$records$metrics$train)
fname.model <- get_model_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range, 
                                   n_layer, n_hidden, n_train)
luz_save(dnn.fit, fname.model)
cat(paste("\n", fname.model, " saved.", sep = ""))
cat("\nSaving model... Done.")


# Evaluate DNN performance w/ predictions 
preds           <- predict(dnn.fit, test_dl)$to(device = "cpu") # get DNN preds
vec.pred.lambda <- as.matrix(preds)[,1] # extract lambda predictions 
vec.pred.mu     <- as.matrix(preds)[,2] # extract mu predictions 
vec.true.lambda <- df$lambda[test_indices] # corresponding lambda true values 
vec.true.mu     <- df$mu[test_indices] # corresponding mu true values 

# Prepare plots 
name.list <- list("lambda", "mu")
true.list <- list("lambda" = vec.true.lambda, "mu" = vec.true.mu)
pred.list <- list("lambda" = vec.pred.lambda, "mu" = vec.pred.mu)

# Save neural network predictions 
save_predictions(pred.list, true.list, nn_type, n_trees, n_taxa,
                  lambda_range, epsilon_range, n_test, n_layer, 
                  n_hidden, n_train)


# Plot predictions (DNN and MLE)
trees_test <- trees[test_indices] # test trees
plot_together_nn_mle_predictions(pred.list, true.list, trees_test, nn_type, n_trees, n_taxa, 
                                 lambda_range, epsilon_range, n_layer, n_hidden, n_train,
                                 save = TRUE)
