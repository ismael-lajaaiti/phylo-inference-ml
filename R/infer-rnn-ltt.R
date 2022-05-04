# Importing libraries 

library(torch)
library(luz)
source("neural-network-functions.R")
source("infer-general-functions.R")


device <- "cuda:1" # GPU where to run computations 
nn_type <- "rnn-ltt" # type of the model: Recurrent Neural Network w/ LTT
model_type <- "bisse"

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
param.range <- param.range.list[[model_type]]

ss_check <- TRUE

# Generate the trees and save 

if (model_type == "crbd"){
  out        <- readPhylogeny(n_trees, n_taxa, param.range, load_trees = TRUE)
  trees      <- out$trees # contains the phylogenetic trees generated
  true.param <- out$param # contains the true values of the parameters 
  
  # Generate LTT data.frame
  df.ltt <- generate_ltt_dataframe(trees, n_taxa, true.param)$ltt
  ds.ltt <- convert_ltt_dataframe_to_dataset(df.ltt, true.param, nn_type)
}

if (model_type == "bisse"){
  out   <- load_dataset_trees(n_trees, n_taxa, param.range, ss_check = ss_check, 
                              load_trees = FALSE)
  true.param <- out$param # contains the true values of the parameters 
  true.param <- true.param[-c(2,3,4,6)] # remove redundant parameters
  
  # Load LTT data.frame
  fname.ltt <- get_dataset_save_name(n_trees, n_taxa, param.range, ss_check)$ltt
  df.ltt <- readRDS(fname.ltt)
}

# Parameters of the NN's training
n_train    <- 90000
n_valid    <- 5000
n_test     <- 5000
n_epochs   <- 100 
batch_size <- 64
patience   <- 10

# Creation of the train, valid and test dataset
train_indices     <- sample(1:n_trees, n_train)
not_train_indices <- setdiff(1:n_trees, train_indices)
valid_indices     <- sample(not_train_indices, n_valid)
test_indices      <- sample(setdiff(not_train_indices, valid_indices), n_test)



if (length(n_taxa) == 1){
  
  # Creation of the dataset 
  train_ds <- ds.ltt(df.ltt[, train_indices], extract_elements(true.param, train_indices))
  valid_ds <- ds.ltt(df.ltt[, valid_indices], extract_elements(true.param, valid_indices))
  test_ds  <- ds.ltt(df.ltt[, test_indices] , extract_elements(true.param, test_indices))
  
  # Creation of the dataloader 
  train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
  valid_dl <- valid_ds %>% dataloader(batch_size=batch_size, shuffle=FALSE)
  test_dl  <- test_ds  %>% dataloader(batch_size=1,          shuffle=FALSE)
}

if (length(n_taxa) == 2){
  
  # Training set 
  true.param.train <- true.param %>% as.data.frame()
  true.param.train <- true.param.train[train_indices, ] %>% as.list()
  list.indices.train <- find_indices_same_size(df.ltt[, train_indices], n_taxa)
  train.set   <- create_all_batch(df.ltt[, train_indices],
                                    true.param.train, list.indices.train,
                                    n_taxa)
  train.set <- reformat_set(train.set, max_batch_size = batch_size)
  
  # Validation set 
  true.param.valid <- true.param %>% as.data.frame()
  true.param.valid <- true.param.valid[valid_indices, ] %>% as.list()
  list.indices.valid <- find_indices_same_size(df.ltt[, valid_indices], n_taxa)
  valid.set   <- create_all_batch(df.ltt[, valid_indices],
                                    true.param.valid, list.indices.valid,
                                    n_taxa)
  valid.set <- reformat_set(valid.set, max_batch_size = batch_size)
  
  # Testing set 
  true.param.test <- true.param %>% as.data.frame()
  true.param.test <- true.param.test[test_indices, ] %>% as.list()
  list.indices.test <- list()
  for (i in 1:n_test){list.indices.test[[i]] = i}
  #list.indices.test <- find_indices_same_size(df.ltt[, test_indices], n_taxa)
  test.set   <- create_all_batch(df.ltt[, test_indices],
                                   true.param.test, list.indices.test,
                                   n_taxa)
  test.set <- reformat_set(test.set, max_batch_size = 1)
}


# Parameters of the RNN
n_hidden  <- 500  # number of neurons in hidden layers 
n_layer   <- 2  # number of stacked RNN layers 
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
opt <- optim_adam(params = rnn$parameters) # optimizer 


# Prepare training 

train_batch <- function(b){
  opt$zero_grad()
  #if (model_type == "crbd"){b$x <- b$x$unsqueeze(2)}
  b$x <- b$x$squeeze(2)$unsqueeze(3)
  output <- rnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$backward()
  opt$step()
  loss$item()
}

valid_batch <- function(b) {
  #if (model_type == "crbd"){b$x <- b$x$unsqueeze(2)}
  b$x <- b$x$squeeze(2)$unsqueeze(3)
  output <- rnn(b$x$to(device = device))
  target <- b$y$to(device = device)
  loss <- nnf_mse_loss(output, target)
  loss$item()
}



# Training loop 

epoch  <- 1
trigger   <- 0 
last_loss <- 100
n_train_batch <- length(train.set)
n_valid_batch <- length(valid.set)

start_time <- Sys.time()

while (epoch <= n_epochs & trigger < patience) {
  
  # Training part 
  rnn$train()
  train_loss <- c()
  
  if (length(n_taxa) == 2){
    random_iter <- sample(1:n_train_batch, n_train_batch)
    c <- 0
    coro::loop(for (i in random_iter) {
      b <- train.set[[i]]
      #print(dim(b$x))
      loss <- train_batch(b)
      train_loss <- c(train_loss, loss)
      c <- c + 1
      #print(c)
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
      b <- valid.set[[i]]
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

end_time <- Sys.time()

run_time <- end_time - start_time
print(run_time)

# Evaluation of the predictions of the RNN w/ test set 

rnn$eval()
nn.pred <- vector(mode = "list", length = n_out)
names(nn.pred) <- names(true.param)


# Compute predictions 

if (length(n_taxa) == 2){
  coro::loop(for (i in 1:n_test) {
    b   <- test.set[[i]]
    out <- rnn(b$x$unsqueeze(3)$to(device = device))
    out <- out$squeeze(1)$to(device = "cpu") %>% as.numeric()
    #true <- b$y %>% as.array()
    for (i in 1:length(out)){
      nn.pred[[i]] <- c(nn.pred[[i]],out[i])
      #test.true[[i]] <- c(test.true[[i]],true[i])
      }
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
fname.mle <- getSaveName(n_trees, n_taxa, param.range)$mle
mle.pred <- readRDS(fname.mle)
mle.pred.test <- as.list(as.data.frame(do.call(cbind, mle.pred))[test_indices,])
if (model_type == "bisse"){mle.pred.test <- mle.pred.test[-c(2,3,4,6)]}
pred.param.test <- list("mle" = mle.pred.test)
pred.param.test[["rnn_ltt"]] <- nn.pred
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
