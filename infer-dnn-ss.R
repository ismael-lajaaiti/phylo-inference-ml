# Importing libraries 

library(torch)
library(luz)
source("summary-statistics.R")
source("infer-general-functions.R")
source("neural-network-functions.R")


nn_type <- "dnn-ss" # type of the model: Deep Neural Network w/ Summary Statistics
model_type <- "bisse" # type of diversification model 

# Parameters of phylogenetic trees
n_trees <- 10000 # total number of trees (train + valid + test)
n_taxa  <- c(100,1000) # size of the trees

# Parameter range of Constant Rate Birth Death model 
param.range.crbd <- list("lambda"  = c(0.1,1.),  # speciation rate
                         "epsilon" = c(0.,0.9))  # extinction rate 

# Parameter range of BiSSE model 
param.range.bisse <- list("lambda"  = c(0.1,1.),  # speciation rate
                          "q"       = c(.01,.1))  # transition rate 0 <-> 1

param.range.list <- list("crbd"  = param.range.crbd,
                         "bisse" = param.range.bisse)

param.range <- param.range.list[[model_type]] # select the range of parameters that 
# correspond to the selected model 
#param.range <- list("lambda" = c(0.1,1.), "epsilon" = c(0.,.9))
ss_check <- TRUE
save_model <- FALSE
save_preds <- FALSE

# Generate the trees and save 
out   <- load_dataset_trees(n_trees, n_taxa, param.range, ss_check = ss_check)
trees      <- out$trees # contains the phylogenetic trees generated 
true.param <- out$param
true.param <- true.param[-c(2,3,4,6)]
name.param <- names(true.param)
n_param    <- length(name.param)

# Create the corresponding summary statistics data.frame
df <- load_dataset_summary_statistics(n_trees, n_taxa, param.range)
df <- df_add_tipstate(df, trees)
df <- scale_summary_statistics(df, n_taxa, name.param)

# Parameters of the NN's training
n_train    <- 9000
n_valid    <- 500
n_test     <- 500
batch_size <- 64

# Creation of the train, valid and test dataset
ds <- convert_ss_dataframe_to_dataset(df)
train_indices <- sample(1:nrow(df), n_train)
not_train_indices <- setdiff(1:nrow(df), train_indices)
valid_indices <- sample(not_train_indices, n_valid)
test_indices  <- setdiff(not_train_indices, valid_indices)
train_ds <- ds(df[train_indices, ], name.param, c("lambda1", "mu0", "mu1", "q10"))
valid_ds <- ds(df[valid_indices, ], name.param, c("lambda1", "mu0", "mu1", "q10"))
test_ds  <- ds(df[test_indices, ], name.param, c("lambda1", "mu0", "mu1", "q10"))

# Creation of the dataloader 
train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
test_dl  <- test_ds  %>% dataloader(batch_size=batch_size, shuffle=FALSE)

# DNN parameters 
n_in      <- length(train_ds[1]$x) # number of neurons of the input layer 
n_out     <- n_param
n_hidden  <- 100 # number of neurons in the hidden layers 
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
    self$fc5 <- nn_linear(in_features = n_hidden, out_features = n_out)
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


if (save_model){
  cat("\nSaving model...")
  n_layer <- length(dnn()$parameters)/2 - 2 
  # epochs <- length(dnn.fit$records$metrics$train)
  fname.model <- get_model_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range, 
                                     n_layer, n_hidden, n_train)
  luz_save(dnn.fit, fname.model)
  cat(paste("\n", fname.model, " saved.", sep = ""))
  cat("\nSaving model... Done.")
}

# Evaluate DNN performance w/ predictions 
preds   <- predict(dnn.fit, test_dl)$to(device = "cpu") # get DNN preds
nn.pred <- preds %>% as.matrix %>% as.data.frame() %>% as.list()
names(nn.pred) <- name.param
if (save_preds){
  # Save neural network predictions 
  save_predictions(pred.list, true.list, nn_type, n_trees, n_taxa,
                    lambda_range, epsilon_range, n_test, n_layer, 
                    n_hidden, n_train)
}


# Prepare plot 
name.param.plot <- c("lambda", "q")

# Plot Predictions 
true.param.test <- as.list(as.data.frame(do.call(cbind, true.param))[test_indices,])
fname.mle <- get_mle_preds_save_name(n_trees, n_taxa, param.range, ss_check)
mle.pred <- readRDS(fname.mle)
mle.pred.test <- as.list(as.data.frame(do.call(cbind, mle.pred))[test_indices,])
mle.pred.test <- mle.pred.test[-c(2,3,4,6)]
pred.param.test <- list("mle" = mle.pred.test)
pred.param.test[[nn_type]] <- nn.pred
param.range.ajusted <- list("lambda" = c(0.1,1.), "q" = c(0.,.1))


plot_pred_vs_true_all(pred.param.test, true.param.test, name.param, param.range.ajusted, 
                      param.range.ajusted, fname = "pred_true_dnnSS_vs_MLE", 
                      save = TRUE)

plot_error_barplot_all(pred.param.test, true.param.test, param.range.ajusted, 
                       save = TRUE, fname = "error_dnnSS_vs_MLE")


