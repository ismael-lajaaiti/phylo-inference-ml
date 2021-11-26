# Importing libraries 

library(torch)
library(luz)
source("neural-network-functions.R")
source("infer-general-functions.R")


device <- "cuda:2" # GPU where to run computations 

nn_type <- "rnn-ltt" # type of the model: Recurrent Neural Network w/ LTT

# Parameters of phylogenetic trees
n_trees <- 1000 # total number of trees (train + valid + test)
n_taxa  <- c(100,1000) # size of the trees
a_range <- c(.1, .5)
b_range <- c(.1, .5)
c_range <- c(.1, .5)
epsilon_range <- c(0, .9)
param.range <- list("a"       = a_range,
                    "b"       = b_range,
                    "c"       = c_range,
                    "epsilon" = epsilon_range)
ss_check <- TRUE

# Generate the trees and save 
out   <- load_dataset_trees(n_trees, n_taxa, param.range, ss_check = ss_check)
trees      <- out$trees # contains the phylogenetic trees generated
true.param <- out$param # contains the true values of the parameters 

# Create the corresponding summary statistics data.frame
out       <- generate_ltt_dataframe(trees, n_taxa, true.param)
df.ltt    <- out$ltt   # ltt dataframe 
df.rates  <- out$rates # rate dataframe
ds.ltt    <- convert_ltt_dataframe_to_dataset(df.ltt, df.rates, nn_type)

# Parameters of the NN's training
n_train    <- 900
n_valid    <- 50
n_test     <- 50
n_epochs   <- 100 
batch_size <- 32
patience   <- 10

# Creation of the train, valid and test dataset
train_indices     <- sample(1:n_trees, n_train)
not_train_indices <- setdiff(1:n_trees, train_indices)
valid_indices     <- sample(not_train_indices, n_valid)
test_indices      <- setdiff(not_train_indices, valid_indices)
train_ds <- ds.ltt(df.ltt[, train_indices], df.rates[train_indices, ])
valid_ds <- ds.ltt(df.ltt[, valid_indices], df.rates[valid_indices, ])
test_ds  <- ds.ltt(df.ltt[, test_indices] , df.rates[test_indices, ])


if (length(n_taxa) == 1){
  # Creation of the dataloader 
  train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
  valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
  test_dl  <- test_ds  %>% dataloader(batch_size=1,          shuffle=FALSE)
}


# Parameters of the RNN
n_hidden  <- 800  # number of neurons in hidden layers 
n_layer   <- 2   # number of stacked RNN layers 
p_dropout <- .01 # dropout probability
n_out     <- length(true.param)


# Build the RNN 

rnn.net <- nn_module(
  initialize = function(n_input, n_out, n_hidden, n_layer, p_dropout = .01,
                        batch_first = TRUE) {
    self$rnn <- nn_lstm(input_size = n_input, hidden_size = n_hidden, 
                        dropout = p_dropout, num_layers = n_layer,
                        batch_first = batch_first)
    self$out <- nn_linear(n_hidden, n_out)
  },
  
  forward = function(x) {
    x <- self$rnn(x)[[1]]
    x <- x[, dim(x)[2], ]
    x %>% self$out() 
  }
)

rnn <- rnn.net(1, n_out, n_hidden, n_layer, p_dropout) # create the RNN
rnn$to(device = device) # move the RNN to the choosen GPU 


find_indices_same_size <- function(df.ltt, n_taxa){
  list.indices <- list()
  for (n in n_taxa[1]:n_taxa[2]){
    list.indices[[as.character(n)]] <- c()
  }
  for(i in 1:ncol(df.ltt)){
      n <- nrow(na.omit(df.ltt[i]))
      list.indices[[as.character(n)]] <- c(list.indices[[as.character(n)]], i)
  }
  return(list.indices)
}


create_one_batch <- function(df.ltt, df.rates, indices){
  
  batch <- list("x" = NA, "y" = NA)
  
  inputs <- na.omit(df.ltt[indices]) %>%
    as.matrix() %>% t() %>%
    torch_tensor()
  
  targets <- df.rates[indices,] %>%
    as.matrix %>%
    torch_tensor()
  
  batch$x <- inputs
  batch$y <- targets
  
  return(batch)
  
}


create_all_batch <- function(df.ltt, df.rates, list.indices, n_taxa){
  
  list.batch <- list()
  #n_iter <- sample(n_taxa[1]:n_taxa[2], length(n_taxa[1]:n_taxa[2]))
  for (n in n_taxa[1]:n_taxa[2]){
    indices <- list.indices[[as.character(n)]]
    batch <- create_one_batch(df.ltt, df.rates, indices)
    if (dim(batch$x)[1] != 0){
      list.batch[[as.character(n)]] <- batch
    }
  }
  return(list.batch)
}



reformat_test_batch <- function(test.batch){
  
  new.test.batch <- list()
  n_test_batch <- length(test.batch)
  
  for (i in 1:n_test_batch){
    batch <- test.batch[[i]]
    n <- dim(batch$x)[1]
    for (j in 1:n){
      mini.batch <- list("x" = NA, "y" = NA)
      input  <- batch$x[j,]
      target <- batch$y[j,]
      mini.batch$x <- input
      mini.batch$y <- target 
      new.test.batch[[length(new.test.batch) + 1]] <- mini.batch
    }
  }
  
  return(new.test.batch)
}


recover_trees_indices <- function(test_indices, list.indices){
  indices <- c()
  for (i in test_indices){
    mapping <- names(list.batch)
    ind <- list.indices[[mapping[i]]]
    indices <- c(indices, ind)
  }
  return(indices)
}


if (length(n_taxa) == 2){
  list.indices <- find_indices_same_size(df.ltt, n_taxa)
  list.batch <- create_all_batch(df.ltt, df.rates, list.indices, n_taxa)
  n_batch <- length(list.batch)
  n_train_batch <- as.integer(n_batch*0.9)
  n_valid_batch <- as.integer(n_batch*0.05)
  n_test_batch  <- n_batch - n_train_batch - n_valid_batch
  train_indices     <- sample(1:n_batch, n_train_batch)
  not_train_indices <- setdiff(1:n_batch, train_indices)
  valid_indices     <- sample(not_train_indices, n_valid_batch)
  test_indices      <- setdiff(not_train_indices, valid_indices)
  recover_test_indices <- recover_trees_indices(test_indices, list.indices)
  train.batch <- list()
  valid.batch <- list()
  test.batch  <- list()
  for (i in 1:n_train_batch){
    ind <- train_indices[i]
    train.batch[[i]] <- list.batch[[ind]]
  }
  for (i in 1:n_valid_batch){
    ind <- valid_indices[i]
    valid.batch[[i]] <- list.batch[[ind]]
  }
  for (i in 1:n_test_batch){
    ind <- test_indices[i]
    test.batch[[i]] <- list.batch[[ind]] 
  } 
  test.batch <- reformat_test_batch(test.batch)
}

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

epoch  <- 1
trigger   <- 0 
last_loss <- 100

while (epoch <= n_epochs & trigger < patience) {
  
  # Training part 
  rnn$train()
  train_loss <- c()
  
  if (length(n_taxa) == 2){
    random_iter <- sample(1:n_train_batch, n_train_batch)
    coro::loop(for (i in random_iter) {
      b <- train.batch[[i]]
      loss <- train_batch(b)
      train_loss <- c(train_loss, loss)
    })
  }
  
  else{
    coro::loop(for (b in train_dl) {
      loss <- train_batch(b)
      train_loss <- c(train_loss, loss)
    })
  }
  
  
  cat(sprintf("epoch %0.3d/%0.3d - train - loss: %3.5f \n",
              epoch, n_epochs, mean(train_loss)))
  
  # Evaluation part 
  rnn$eval()
  valid_loss <- c()
  
  
  if (length(n_taxa) == 2){
    coro::loop(for (i in 1:n_valid_batch) {
      b <- valid.batch[[i]]
      loss <- valid_batch(b)
      valid_loss <- c(valid_loss, loss)
    })
  }
  
  else{
    coro::loop(for (b in valid_dl) {
      loss <- valid_batch(b)
      valid_loss <- c(valid_loss, loss)
    })
  }
  
  
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
cat("\nSaving model...")
fname.model <- get_model_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range, 
                                   n_layer, n_hidden, n_train)
torch_save(rnn, fname.model)
cat(paste("\n", fname.model, " saved.", sep = ""))
cat("\nSaving model... Done.")

# Evaluation of the predictions of the RNN w/ test set 

rnn$eval()
nn.pred <- vector(mode = "list", length = n_out)
names(nn.pred) <- names(true.param)

# Compute predictions 

if (length(n_taxa) == 2){
  random_iter <- sample(1:length(test.batch), 40)
  test_indices <- recover_test_indices[random_iter]
  coro::loop(for (i in random_iter) {
    b <- test.batch[[i]]
    out <- rnn(b$x$unsqueeze(2)$unsqueeze(1)$to(device = device))
    out <- out$squeeze(1)$to(device = "cpu") %>% as.numeric()
    for (i in 1:length(out)){nn.pred[[i]] <- c(nn.pred[[i]],out[i])}
  })
}

if (length(n_taxa) == 1) {
  coro::loop(for (b in test_dl) {
    out <- rnn(b$x$reshape(c(b$x$shape, 1L))$to(device = device))
    pred <- as.numeric(out$to(device = "cpu")) # move the tensor to CPU 
    true <- as.numeric(b$y)
    vec.pred.lambda <- c(vec.pred.lambda, pred[1])
    vec.pred.mu     <- c(vec.pred.mu, pred[2])
    vec.true.lambda <- c(vec.true.lambda, true[1])
    vec.true.mu     <- c(vec.true.mu, true[2])
  })
}


# Prepare plot 
true.param.test <- as.list(as.data.frame(do.call(cbind, true.param))[test_indices,])
fname.mle <- get_mle_preds_save_name(n_trees, n_taxa, param.range, ss_check)
mle.pred <- readRDS(fname.mle)
mle.pred.test <- as.list(as.data.frame(do.call(cbind, mle.pred))[test_indices,])
pred.param.test <- list("mle" = mle.pred.test)
pred.param.test[[nn_type]] <- nn.pred
param.range.ajusted <- param.range[-4]
param.range.ajusted[["mu"]] <- c(param.range[["c"]][1]*param.range[["epsilon"]][1],
                                 param.range[["c"]][2]*param.range[["epsilon"]][2])

plot_pred_vs_true_all(pred.param.test, true.param.test, name.param, param.range.ajusted)


# Plot Predictions 
trees_test <- trees[test_indices] # test trees
plot_together_nn_mle_predictions(pred.list, true.list, trees_test, nn_type, n_trees, n_taxa, 
                                 lambda_range, epsilon_range, n_layer, n_hidden, n_train,
                                 save = FALSE)


pred.list.mle <- get_mle_preds(trees_test)
pred.list.all <- list()
pred.list.all[[nn_type]] <- pred.list
pred.list.all[["mle"]]   <- pred.list.mle

plot_bars_mle_vs_nn(pred.list.all, true.list, nn_type, name.list, save = FALSE, n_trees, n_taxa, 
                    lambda_range, epsilon_range, n_test, n_layer, n_hidden, n_train)
