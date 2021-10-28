# This file contains the general functions needed during the interference 
# processes 


#### Libraries & Sources ####

library(ape)
library(diversitree)
library(RPANDA)
library(MLmetrics)
source("summary-statistics.R")

#### end ####



#' Compute r and epsilon
#'
#' Compute diversification rate (r) and turnover rate (epsilon) given the 
#' speciation rate (lambda) and extinction rate (lambda). 
#' Formally : 
#' r = lambda - mu
#' epsilon = mu / lambda
#'
#' @param vec.lambda_mu vector of length 2 containing lambda and mu
#'
#' @return vec.r_epsilon vector of length 2 containing r and epsilon 
#' @export
#' @examples
get_r_epsilon_single <- function(vec.lambda_mu){
  lambda  <- vec.lambda_mu[1]
  mu      <- vec.lambda_mu[2]
  r       <- lambda - mu
  epsilon <- mu / lambda
  vec.r_epsilon <- c(r, epsilon)
  return(vec.r_epsilon)
}


#' Compute lambda and mu 
#'
#' Compute speciation rate (lambda) and extinction rate (mu) given the 
#' diversification rate (r) and turnover rate (epsilon). 
#' Formally : 
#' lambda  =            r / (1 - epsilon)
#' mu      = epsilon *  r / (1 - epsilon)
#'
#' @param vec.r_epsilon vector of length 2 containing r and epsilon
#'
#' @return vec.lambda_mu vector of length 2 containing lambda and mu 
#' @export
#' @examples
get_lambda_mu_single <- function(vec.r_epsilon){
  r       <- vec.r_epsilon[1]
  epsilon <- vec.r_epsilon[2]
  lambda  <- r / (1 - epsilon)
  mu      <- epsilon * lambda 
  vec.lambda_mu <- c(lambda, mu)
  return(vec.lambda_mu)
}


#' Extend get_r_epsilon_single to vectors  
#'
#' Given a vector of speciation rates and a vector of extinction rates
#' compute the corresponding diversification rates and turnover rates
#'
#' @param list.lambda_mu $lambda vector of speciation rates 
#'                       $mu     vector of extinction rates
#'
#' @return list.r_epsilon $r        vector of diversification rates
#'                        $epsilon  vector of turnover rates 
#' @export
#' @examples
get_r_epsilon_extended <- function(list.lambda_mu){
  vec.lambda  <- list.lambda_mu$lambda
  vec.mu      <- list.lambda_mu$mu
  vec.r       <- c()
  vec.epsilon <- c()
  n           <- length(vec.mu)
  
  for (i in 1:n){
    vec.lambda_mu <- c(vec.lambda[i], vec.mu[i])
    vec.r_epsilon <- get_r_epsilon_single(vec.lambda_mu)
    r <- vec.r_epsilon[1]
    epsilon <- vec.r_epsilon[2]
    vec.r <- c(vec.r, r)
    vec.epsilon <- c(vec.epsilon, epsilon)
  }
  
  list.r_epsilon <- list("r" = vec.r, "epsilon" = vec.epsilon)
  
  return(list.r_epsilon)
}


#' Extend get_lambda_mu_single to vectors  
#'
#' Given a vector of diversification rates and a vector of turnover rates
#' compute the corresponding speciation rates and extinction rates
#'
#' @param list.r_epsilon $r        vector of diversification rates
#'                        $epsilon  vector of turnover rates 
#'
#' @return list.lambda_mu $lambda vector of speciation rates 
#'                       $mu     vector of extinction rates
#' @export
#' @examples
get_lambda_mu_list <- function(list.r_epsilon){
  vec.r       <- list.r_epsilon$r
  vec.epsilon <- list.r_epsilon$epsilon
  vec.lambda  <- c()
  vec.mu      <- c()
  n           <- length(vec.r)
  
  for (i in 1:n){
    vec.r_epsilon <- c(vec.r[i], vec.epsilon[i])
    vec.lambda_mu <- get_lambda_mu_single(vec.r_epsilon)
    lambda <- vec.lambda_mu[1]
    mu     <- vec.lambda_mu[2]
    vec.lambda <- c(vec.lambda, lambda)
    vec.mu     <- c(vec.mu, mu)
  }
  
  list.lambda_mu <- list("lambda" = vec.lambda, "mu" = vec.mu)
  
  return(list.lambda_mu)
}


get_mu_vec <- function(list.lambda_r){
  vec.lambda <- list.lambda_r$lambda
  vec.r      <- list.lambda_r$r
  vec.mu     <- c()
  n          <- length(vec.r)
  
  for (i in 1:n){
    mu <- vec.lambda[i] - vec.r[i]
    vec.mu <- c(vec.mu, mu)
  }
  
  return(vec.mu)
  
}



plot_pred_vs_true <- function(pred.list, true.list, names,
                              method = "model", save = FALSE, fname = "file-name",
                              r2_score = TRUE, lm_fit = TRUE, subtitles = ""){
  
  n <- length(pred.list)
  n_row <- n %/% 2
  n_col <- 2
  
  if (save){
    fname <- paste(fname, "pdf", sep=".")
    aspect_ratio <- 1.62 * (2 / n_row)
    pdf(fname, width = 10, height = 10/aspect_ratio, pointsize = 15/sqrt(aspect_ratio))
    }
  
  par(mfrow=c(n_row,n_col))
  
  for (i in 1:n){
    
    pred.name <- paste(names[[1 + (i-1)%%2]], method, sep=".") # parameter name
    pred <- pred.list[[i]] # parameter prediction 
    true <- true.list[[1 + (i-1)%%2]] # parameter true 
    
    # Evaluate the R2 score of predictions vs. truth
    if (r2_score){
      r2 <- R2_Score(pred, true) # compute r2
      r2 <- format(round(r2, 3), nsmall = 3) # format r2 
      plot(true, pred, main = paste(names[[1 + (i-1)%%2]], "- r2 =", r2, sep=" "),
           sub = subtitles[1 + (i-1)%/%2])
    }
    else{
      plot(true, pred, main = names[[i]])
    }
    
    abline(0, 1) # plot identity line (for the eye)
    
    # Linear fit of predictions vs. truth (to see significant trends)
    if (lm_fit){
      fit = lm(pred ~ true)
      sig = summary(fit)$coefficients[2,4]
      abline(fit, col="red", lty = ifelse(sig < .05,1,2))
    }
  }
  
  if (save){dev.off()}
}


generate_trees <- function(n_trees, n_taxa, lambda_range, epsilon_range,
                           ss_check = TRUE){
  
  trees <- list() # initialize tree list where trees will be stored 
  vec.true.lambda <- c() # vector where true lambda values will be stored
  vec.true.mu     <- c() # vector where true mu values will be stored 
  
  cat("Generation of trees...\n")
  
  while (length(trees) < n_trees){
    # Generate the phylogenetic tree 
    true.lambda  <- runif(1, lambda_range[1] , lambda_range[2]) # generate random lambda
    epsilon      <- runif(1, epsilon_range[1], epsilon_range[2]) # generate random epsilon
    true.mu      <- epsilon * true.lambda # compute corresponding mu
    tree <- trees(c(true.lambda, true.mu), "bd", max.taxa=n_taxa)[[1]] # create the tree 
    
    # If checking that summary statistics have no NA
    if (ss_check){
      ss <- get_ss(tree) # compute summary statistics 
      if (!any(is.na(ss))){ # if no NA values 
        trees <- append(trees, list(tree)) # append the new tree to the list
        vec.true.lambda  <- c(vec.true.lambda, true.lambda) # store lambda prediction
        vec.true.mu      <- c(vec.true.mu, true.mu) # store mu prediction 
        progress(length(trees), n_trees, progress.bar = TRUE,
                 init = (length(trees)==1)) # print progression
        }
    }
    
    # Else just append the new tree to the list and save rtes 
    else{
      trees <- append(trees, list(tree))
      vec.true.lambda  <- c(vec.true.lambda, true.lambda) # store lambda prediction
      vec.true.mu      <- c(vec.true.mu, true.mu) # store mu prediction 
      progress(length(trees), n_trees, progress.bar = TRUE,
               init = (length(trees)==1)) # print progression
      }
    
  }
  
  cat("\nGeneration of trees... Done.")
  
  # Prepare output containing: trees list, true lambda vector, true mu vector.
  out <- list("trees"  = trees, 
              "lambda" = vec.true.lambda,
              "mu"     = vec.true.mu)
  
  return(out)
  
}


generate_ltt_dataframe <- function(trees, n_taxa, vec.lambda, vec.mu){
  
  n_trees  <- length(trees) # number of trees 
  df.ltt   <- data.frame("tree1" = rep(NA,n_taxa)) # initialize df.ltt
  df.rates <- data.frame("lambda" = rep(NA, n_trees), # initialize df.rates
                         "mu" = rep(NA, n_trees)) 
  df.rates$lambda <- vec.lambda # fill lambda column w/ corresponding values
  df.rates$mu     <- vec.mu # fill mu column w/ corresponding values 
  
  cat("Creating LTT dataframe...\n")
  
  for (i in 1:n_trees){
    tree <- trees[[i]] # get tree 
    ltt.coord <- ape::ltt.plot.coords(tree) # get ltt coordinates 
    ltt.coord <- as.data.frame(ltt.coord)
    df.ltt[paste("tree", i, sep = "")] <- ltt.coord$time
    progress(i, n_trees, progress.bar = TRUE, init = (i==1))
  }
  
  cat("\nCreating LTT dataframe... Done.")

  out <- list("ltt" = df.ltt, "rates" = df.rates) # function output
  
  return(out)
  
}

convert_ltt_dataframe_to_dataset <- function(df.ltt, df.rates){
  
  ds.ltt <- torch::dataset(
    name <- "ltt_dataset", 
    initialize = function(df.ltt, df.rates){
      
      array.ltt <- df.ltt %>% 
        as.matrix() %>% 
        array(dim = c(n_taxa, ncol(df.ltt), 1)) # convert df to array of good dimension
      
      # input data 
      x <- array.ltt
      self$x <- x
      
      # target data 
      target.names <- c("lambda", "mu")
      y = df.rates[target.names] %>% 
        as.matrix()
      self$y <- torch_tensor(y)
      
    }, 
    
    .getitem = function(i) {
      list(x = self$x[i, , ,drop=TRUE], y = self$y[i, ])
    }, 
    
    .length = function() {
      self$y$size()[[1]]
    }
  )
  
  return(ds.ltt)
}


plot_mle_predictions <- function(trees, vec.true.lambda, vec.true.mu, 
                                 lambda_range, mu_range, type = "all", 
                                 save = FALSE, fname = NA){
  
  n_trees <- length(trees)
  n_taxa  <- length(trees[[1]]$tip.label)
  
  vec.pred.rpanda.lambda      <- c()
  vec.pred.rpanda.mu          <- c()
  vec.pred.diversitree.lambda <- c()
  vec.pred.diversitree.mu     <- c()
  vec.pred.ape.lambda         <- c()
  vec.pred.ape.mu             <- c()
  
  
  cat("Compute MLE predictions...\n")
  
  for (i in 1:n_trees){
    
    tree <- trees[[i]] # get the phylogenetic tree 
    
    # Prepare the fit 
    height <- max(node.age(tree)$ages) # height of the tree 
    f.lambda <-function(t,y){y[1]} # lambda is a function constant over time 
    f.mu     <-function(t,y){y[1]} # mu is a function constant over time 
    lambda_par <- c(.2)  # initial value for lambda
    mu_par     <- c(.1)  # initial value for mu
    
    # RPANDA - Infer lambda & mu 
    if (type == "all" | type == "rpanda"){
      fit.rpanda <- RPANDA::fit_bd(tree, height, f.lambda, f.mu, lambda_par, mu_par, 
                                   cst.lamb = TRUE, cst.mu = TRUE) # fit bd model 
      pred.rpanda.lambda <- fit.rpanda$lamb_par # get lambda prediction 
      pred.rpanda.mu     <- fit.rpanda$mu_par # get mu prediction 
      vec.pred.rpanda.lambda <- c(vec.pred.rpanda.lambda,
                                  pred.rpanda.lambda) # store lambda prediction
      vec.pred.rpanda.mu     <- c(vec.pred.rpanda.mu,
                                  pred.rpanda.mu) # store mu prediction 
    }
    
    # diversitree - Infer lambda & mu
    if (type == "all" | type == "diversitree"){
      lik  <- diversitree::make.bd(tree)
      fit.diversitree <- diversitree::find.mle(lik, c(lambda_par, mu_par)) # fit bd model 
      pred.diversitree.lambda <- fit.diversitree$par[[1]] # get lambda prediction 
      pred.diversitree.mu     <- fit.diversitree$par[[2]] # get mu prediction 
      vec.pred.diversitree.lambda <- c(vec.pred.diversitree.lambda,
                                       pred.diversitree.lambda) # store lambda prediction
      vec.pred.diversitree.mu     <- c(vec.pred.diversitree.mu,
                                       pred.diversitree.mu) # store mu prediction
    }
    
    # ape - Infer r & epsilon then compute lambda & mu 
    if (type == "all" | type == "ape"){
      fit.ape <- ape::birthdeath(tree) # fit bd model 
      pred.ape.r <- fit.ape$para[[2]] # get r prediction 
      pred.ape.epsilon <- fit.ape$para[[1]] # get epsilon prediction
      # compute lambda and mu from r and epsilon 
      pred.ape.lambda_mu <- get_lambda_mu_single(c(pred.ape.r, pred.ape.epsilon))
      pred.ape.lambda <- pred.ape.lambda_mu[1] # get lambda value 
      pred.ape.mu     <- pred.ape.lambda_mu[2] # get mu value 
      vec.pred.ape.lambda <- c(vec.pred.ape.lambda, pred.ape.lambda) # store lambda prediction
      vec.pred.ape.mu <- c(vec.pred.ape.mu, pred.ape.mu) # store mu prediction 
    }
    
    progress(i, n_trees, progress.bar = TRUE, init = (i==1)) # print progression
    
  }
  
  cat("\nCompute MLE predictions... Done.")
  
  names <- list("lambda", "mu")
  true.list <- list(vec.true.lambda, vec.true.mu)
  
  
  if (type == "all"){subtitles <- list("rpanda", "diversitree", "ape")}
  else{subtitles <- list("")}
  
  if (type == "rpanda"){
    pred.list <- list(vec.pred.rpanda.lambda, vec.pred.rpanda.mu)
  }
  
  else if (type == "diversitree"){
    pred.list <- list(vec.pred.diversitree.lambda, vec.pred.diversitree.mu)
  }
  
  else if (type == "ape"){
    pred.list <- list(vec.pred.ape.lambda, vec.diversitree.ape.mu)
  }
  
  else if (type == "all"){
    pred.list <- list(vec.pred.rpanda.lambda, vec.pred.rpanda.mu,
                      vec.pred.diversitree.lambda, vec.pred.diversitree.mu, 
                      vec.pred.ape.lambda, vec.pred.ape.mu)
  }
  
  else{print("Error: type unkown. Type should be either: 'rpanda', 'diversitree', 'ape' or 'all'.")}
  
  
  if (is.na(fname)){
    path <- "figures/mle/"
    fname <- paste("mle", type, "ntaxa", n_taxa,
                   "lambda", lambda_range[1], lambda_range[2], 
                   "mu", mu_range[1], mu_range[2],
                   "ntest", n_trees, sep="-")
    fname <- paste(path, fname, sep="")
  }
  
  plot_pred_vs_true(pred.list, true.list, names, "mle", save = save, fname = fname, 
                    subtitle = subtitles)
  
}


get_dataset_save_names <- function(n_trees, n_taxa, lambda_range, mu_range,
                                   ss_check){
  
  dir <- "trees-dataset/"
  fname <- paste("ntrees", n_trees, "ntaxa", n_taxa,
                 "lambda" , lambda_range[1] , lambda_range[2], 
                 "epsilon", epsilon_range[1], epsilon_range[2], 
                 "sscheck", ss_check, sep="-")
  fname.trees  <- paste(dir, fname, "-trees.rds", sep="")
  fname.lambda <- paste(dir, fname, "-lambda.rds", sep="")
  fname.mu     <- paste(dir, fname, "-mu.rds", sep="")
  
  fnames <- list("trees"  = fname.trees, 
                 "lambda" = fname.lambda, 
                 "mu"     = fname.mu)
  
  return(fnames)
  
}




