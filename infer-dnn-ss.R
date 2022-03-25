# Importing libraries 

library(torch)
library(luz)
source("summary-statistics.R")
source("infer-general-functions.R")
source("neural-network-functions.R")


nn_type <- "dnn-ss" # type of the model: Deep Neural Network w/ Summary Statistics
model_type <- "crbd" # type of diversification model 

# Parameters of phylogenetic trees
n_trees <- 100000 # total number of trees (train + valid + test)
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
out   <- load_dataset_trees(n_trees, n_taxa, param.range, load_trees = TRUE, ss_check = ss_check)
trees      <- out$trees # contains the phylogenetic trees generated 
true.param <- out$param
if (model_type == "bisse"){true.param <- true.param[-c(2,3,4,6)]}
name.param <- names(true.param)
n_param    <- length(name.param)

# Create the corresponding summary statistics data.frame
df <- readSummaryStatistics(n_trees, n_taxa, param.range)
if (model_type == "bisse"){df <- df_add_tipstate(df, trees)}
df <- scale_summary_statistics(df, n_taxa, c("lambda", "mu"))

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
train_ds <- ds(df[train_indices, ], c("lambda", "mu"), c())
valid_ds <- ds(df[valid_indices, ], c("lambda", "mu"), c())
test_ds  <- ds(df[test_indices, ],  c("lambda", "mu"), c())

# Creation of the dataloader 
train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
test_dl  <- test_ds  %>% dataloader(batch_size=1, shuffle=FALSE)

# DNN parameters 
n_in      <- length(train_ds[1]$x) # number of neurons of the input layer 
n_out     <- length(name.param)
n_hidden  <- 100 # number of neurons in the hidden layers 
p_dropout <- 0.01 # dropout probability 
n_epochs  <- 100 # maximum number of epochs for the training 
patience  <- 4 # patience of the early stopping 

# Build the DNN 
dnn.net <- nn_module(
  
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

dnn <- dnn.net() # create CNN
dnn$to(device = device) # Move it to the choosen GPU
opt <- optim_adam(params = dnn$parameters) # optimizer 

train_batch <- function(b){
  opt$zero_grad()
  output <- dnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$backward()
  opt$step()
  loss$item()
}

valid_batch <- function(b) {
  output <- dnn(b$x$to(device = device))
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
  
  # Training 
  dnn$train()
  train_loss <- c()
  coro::loop(for (b in train_dl) { # loop over batches 
    loss <- train_batch(b)
    train_loss <- c(train_loss, loss)
  })
  
  # Print Epoch and value of Loss function 
  cat(sprintf("epoch %0.3d/%0.3d - train - loss: %3.5f \n",
              epoch, n_epochs, mean(train_loss)))
  
  # Validation 
  dnn$eval()
  valid_loss <- c()
  coro::loop(for (b in test_dl) { # loop over batches 
    loss <- valid_batch(b)
    valid_loss <- c(valid_loss, loss)
  })
  current_loss <- mean(valid_loss)
  
  # Early Stopping 
  if (current_loss > last_loss){trigger <- trigger + 1} 
  else{
    trigger   <- 0
    last_loss <- current_loss
  }
  
  # Print Epoch and value of Loss function
  cat(sprintf("epoch %0.3d/%0.3d - valid - loss: %3.5f \n",
              epoch, n_epochs, current_loss))
  
  epoch <- epoch + 1 
}




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



dnn$eval()
nn.pred <- vector(mode = "list", length = n_out)
names(nn.pred) <- names(true.param)

# Compute predictions 
coro::loop(for (b in test_dl) {
  out <- dnn(b$x$to(device = device))
  pred <- as.numeric(out$to(device = "cpu")) # move the tensor to CPU 
  #true <- as.numeric(b$y)
  for (i in 1:n_out){nn.pred[[i]] <- c(nn.pred[[i]], pred[i])}
})


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
name.param.plot <- c("lambda", "mu")

# Prepare plot 
true.param.test <- as.list(as.data.frame(do.call(cbind, true.param))[test_indices,])
fname.mle <- getSaveName(n_trees, n_taxa, param.range)$mle
mle.pred <- readRDS(fname.mle)
mle.pred.test <- as.list(as.data.frame(do.call(cbind, mle.pred))[test_indices,])
if (model_type == "bisse"){mle.pred.test <- mle.pred.test[-c(2,3,4,6)]}
pred.param.test <- list("mle" = mle.pred.test)
pred.param.test[["dnn-ss"]] <- nn.pred
if (model_type == "bisse"){
  param.range.ajusted <- list("lambda" = c(0.1,1.), "q" = c(0.,.1))
  param.range.in      <- list("lambda" = c(0.2,.9), "q" = c(.02,.09))
} else{
  param.range.ajusted <- list("lambda" = c(0.1,1.), "q" = c(0.,1.))
  param.range.in      <- list("lambda" = c(0.2,.9), "q" = c(0.,.8))
}

plot_pred_vs_true_all(pred.param.test, true.param.test, name.param, param.range.ajusted, 
                      param.range.in, fname = "", 
                      save = FALSE)

#reord_names <- c("mle","cnn_cblv", "dnn-ss", "cnn-ltt", "rnn-ltt", "gnn-phylo")
plot_error_barplot_all(pred.param.test, true.param.test, param.range.in, name.param,
                       save = FALSE, fname = "")


