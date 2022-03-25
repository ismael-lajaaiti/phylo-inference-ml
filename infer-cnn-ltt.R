# Importing Libraries and Sources

library(torch)
source("infer-general-functions.R")
source("neural-network-functions.R")


device <- "cuda:1" # GPU where to run computations 

nn_type <- "cnn-ltt" # type of the network: Convolutional Neural Network w/ LTT
model_type <- "bisse" # type of diversification model 

# Parameters of phylogenetic trees
n_trees <- 100000 # total number of trees (train + valid + test)
n_taxa  <- c(100, 1000) # size of the trees


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

ss_check <- TRUE

# Generate the trees and save 
out   <- readPhylogeny(n_trees, n_taxa, param.range, ss_check = ss_check, 
                            load_trees = FALSE)
true.param <- out$param # contains the true values of the parameters 
if(model_type == "bisse"){true.param <- true.param[-c(2,3,4,6)]}

# Create the corresponding summary statistics data.frame
fname.ltt <- getSaveName(n_trees, n_taxa, param.range, ss_check)$ltt
df.ltt <- readRDS(fname.ltt)
df.ltt <- generate_ltt_dataframe(trees, n_taxa, true.param)$ltt
ds.ltt <- convert_ltt_dataframe_to_dataset(df.ltt, true.param, nn_type)

# Parameters of the NN's training
n_train    <- 90000 
n_valid    <- 5000
n_test     <- 5000
n_epochs   <- 100
batch_size <- 64
patience   <- 3

# Creation of the train, valid and test dataset
train_indices     <- sample(1:n_trees, n_train)
not_train_indices <- setdiff(1:n_trees, train_indices)
valid_indices     <- sample(not_train_indices, n_valid)
test_indices      <- setdiff(not_train_indices, valid_indices)

train_ds <- ds.ltt(df.ltt[, train_indices], extract_elements(true.param, train_indices))
valid_ds <- ds.ltt(df.ltt[, valid_indices], extract_elements(true.param, valid_indices))
test_ds  <- ds.ltt(df.ltt[, test_indices] , extract_elements(true.param, test_indices))

# Creation of the dataloader 
train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
test_dl  <- test_ds  %>% dataloader(batch_size=1,          shuffle=FALSE)


n_hidden  <- 8
n_layer   <- 3
ker_size  <- 5
p_dropout <- 0.01
n_input   <- max(n_taxa)
n_out     <- length(true.param)

# Build the CNN

cnn.net <- nn_module(
  
  "corr-cnn",
  
  initialize = function(n_input, n_out, n_hidden, n_layer, ker_size, p_dropout) {
    self$conv1 <- nn_conv1d(in_channels = 1, out_channels = n_hidden, kernel_size = ker_size)
    self$conv2 <- nn_conv1d(in_channels = n_hidden, out_channels = 2*n_hidden, kernel_size = ker_size)
    self$conv3 <- nn_conv1d(in_channels = 2*n_hidden, out_channels = 2*2*n_hidden, kernel_size = ker_size)
    n_flatten <- compute_dim_ouput_flatten_cnn(n_input, n_layer, ker_size)
    self$fc1 <- nn_linear(in_features = n_flatten * (2*2*n_hidden), out_features = 100)
    self$fc2 <- nn_linear(in_features = 100, out_features = n_out)
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
      nnf_relu() %>%
      nnf_dropout(p = p_dropout) %>%
      
      self$fc2()
  }
)

cnn <- cnn.net(n_input, n_out, n_hidden, n_layer, ker_size, p_dropout) # create CNN
cnn$to(device = device) # Move it to the choosen GPU

# Prepare training 

opt <- optim_adam(params = cnn$parameters) # optimizer 

train_batch <- function(b){
  opt$zero_grad()
  #if (model_type == "crbd"){b$x <- b$x$unsqueeze(2)}
  output <- cnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$backward()
  opt$step()
  loss$item()
}

valid_batch <- function(b) {
  #if (model_type == "crbd"){b$x <- b$x$unsqueeze(2)}
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
  
  # Training 
  cnn$train()
  train_loss <- c()
  coro::loop(for (b in train_dl) { # loop over batches 
    loss <- train_batch(b)
    train_loss <- c(train_loss, loss)
  })
  
  # Print Epoch and value of Loss function 
  cat(sprintf("epoch %0.3d/%0.3d - train - loss: %3.5f \n",
              epoch, n_epochs, mean(train_loss)))
  
  # Validation 
  cnn$eval()
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
  #if (model_type == "crbd"){b$x <- b$x$unsqueeze(2)}
  out <- cnn(b$x$to(device = device))
  pred <- as.numeric(out$to(device = "cpu")) # move the tensor to CPU 
  #true <- as.numeric(b$y)
  for (i in 1:n_out){nn.pred[[i]] <- c(nn.pred[[i]], pred[i])}
})

# Prepare plot 
name.param <- ifelse(model_type == "crbd", c("lambda", "mu"), c("lambda", "q"))

# Save neural network predictions 
#save_predictions(pred.list, true.list, nn_type, n_trees, n_taxa,
#                 lambda_range, epsilon_range, n_test, n_layer, 
#                 n_hidden, n_train, ker_size)                                              


# Prepare plot 
true.param.test <- as.list(as.data.frame(do.call(cbind, true.param))[test_indices,])
fname.mle <- getSaveName(n_trees, n_taxa, param.range)$mle
mle.pred <- readRDS(fname.mle)
mle.pred.test <- as.list(as.data.frame(do.call(cbind, mle.pred))[test_indices,])
if (model_type == "bisse"){mle.pred.test <- mle.pred.test[-c(2,3,4,6)]}
pred.param.test <- list("mle" = mle.pred.test)
pred.param.test[["cnn_ltt"]] <- nn.pred
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

ind.list <- list("train" = train_indices,
                 "valid" = valid_indices,
                 "test"  = test_indices)

saveRDS(ind.list, "indices_sets_crbd.rds")


pred_list <- readRDS("predGNN_crbd.rds")
pred <- c()
for (i in 1:5000){
  pred <- c(pred, pred_list[[2]][[i]])
}
predGNN <- list()
predGNN$mu <- pred
par(mfrow = c(1,1))
plot(pred, true.param.test$mu)


names(pred.param.test)
pred.param.test[["gnn_phylo"]] <- predGNN

param.test.save <- list("true" = true.param.test, "pred" = pred.param.test)
saveRDS(param.test.save, "param_bisse.rds")
