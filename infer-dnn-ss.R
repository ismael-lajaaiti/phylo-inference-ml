library(torch)
library(luz)
library(ggplot2)
library(MLmetrics)
source("summary-statistics.R")
source("infer-general-functions.R")
source("neural-network-functions.R")

# Parameters of phylogenetic trees
n_trees <- 5200 # total number of trees (train + test)
n_taxa <- 100 # size of the trees
lambda_range <- c(0.15, 0.25) # range within random lambda will be generated 
mu_range <- c(0.0, 0.1) # same for mu

# Generate the data.frame
df <- generate_ss_dataframe(n_trees, n_taxa,
                            lambda_range[1], lambda_range[2],
                            mu_range[1], mu_range[2])

# Parameters of the NN's training
n_train <- 5000
batch_size <- 64
n_epochs <- 100

# Creation of the test and train dataset
ds <- convert_dataframe_to_dataset(df)
train_indices <- sample(1:nrow(df), n_train) 
test_indices <- setdiff(1:nrow(df), train_indices) 
train_ds <- ds(df[train_indices, ], direct_target=FALSE)
test_ds  <- ds(df[test_indices, ] , direct_target=FALSE)

# Creation of the dataloader 
train_dl <- train_ds %>% dataloader(batch_size=batch_size, shuffle=TRUE)
test_dl  <- test_ds  %>% dataloader(batch_size=batch_size, shuffle=FALSE)

n_in <- ncol(df) - 4
n_hidden <- n_in
p_dropout <- 0.01

net <- nn_module(
  
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

net.fit.1 <- net %>%
  setup(
    loss = function(y_hat, y_true) nnf_mse_loss(y_hat, y_true),
    optimizer = optim_adam
  ) %>%
  fit(train_dl, epochs = n_epochs, valid_data = test_dl, 
      callbacks = list(luz_callback_early_stopping(patience = 7)))

preds <- predict(net.fit.1, test_dl)
preds.r <- as.matrix(preds)[,1]
preds.epsilon <- as.matrix(preds)[,2]
test_dl <- dataloader(test_ds, batch_size = n_trees - n_train)
targets <- (test_dl %>% dataloader_make_iter() %>% dataloader_next())$y %>% 
  as.matrix()
targets.r <- targets[,1]
targets.epsilon <- targets[,2]

train_ds2 <- ds(df[train_indices, ], direct_target=TRUE)
test_ds2  <- ds(df[test_indices, ] , direct_target=TRUE)

train_dl2 <- train_ds2 %>% dataloader(batch_size=batch_size, shuffle=TRUE)
test_dl2  <- test_ds2  %>% dataloader(batch_size=batch_size, shuffle=FALSE)

net.fit.2 <- net %>%
  setup(
    loss = function(y_hat, y_true) nnf_mse_loss(y_hat, y_true),
    optimizer = optim_adam
  ) %>%
  fit(train_dl2, epochs = n_epochs, valid_data = test_dl2, 
      callbacks = list(luz_callback_early_stopping(patience = 7)))

preds <- predict(net.fit.2, test_dl2)
preds.lambda <- as.matrix(preds)[,1]
preds.mu <- as.matrix(preds)[,2]
test_dl2 <- dataloader(test_ds2, batch_size = n_trees - n_train)
targets <- (test_dl2 %>% dataloader_make_iter() %>% dataloader_next())$y %>% 
  as.matrix()
targets.lambda <- targets[,1]
targets.mu <- targets[,2]


list.r_epsilon.targets <- list("r" = targets.r, "epsilon" = targets.epsilon)
list.r_epsilon.preds <- list("r" = preds.r, "epsilon" = preds.epsilon)
list.lambda_mu.targets <- get_lambda_mu_list(list.r_epsilon.targets)
list.lambda_mu.preds <- get_lambda_mu_list(list.r_epsilon.preds)

targets.list <- list(targets.lambda, targets.mu, 
                     targets.r, targets.epsilon, 
                     list.lambda_mu.targets$lambda, list.lambda_mu.targets$mu)
preds.list <- list(preds.lambda, preds.mu, preds.r, preds.epsilon, 
                   list.lambda_mu.preds$lambda, list.lambda_mu.preds$mu)
names <- list("lambda", "mu", "r", "epsilon", "lambda(r,epsilon)", "mu(r,epsilon)")

fname <- paste("taxa", n_taxa, 
               "lambda", lambda_range[1], lambda_range[2],
               "mu", mu_range[1], mu_range[2], 
               "train", n_train, 
               "test", n_trees-n_train, 
               "slopes", 10, sep="-")


path <- "figures/dnn-ss/"

ftype <- ".pdf"

full_fname <- paste(path, fname, ftype, sep="")

pdf(full_fname)



par(mfrow=c(3,2))

for (i in 1:6){
  pred <- preds.list[[i]]
  true <- targets.list[[i]]
  r2.lambda <- R2_Score(pred, true)
  r2.lambda <- format(round(r2.lambda, 3), nsmall = 3)
  plot(true, pred, main=paste(names[[i]], "- r2 =", r2.lambda, sep=" "))
  abline(0,1)
  fit = lm(pred ~ true)
  sig = summary(fit)$coefficients[2,4]
  abline(fit, col="red", lty = ifelse(sig < .05,1,2))
}

dev.off()