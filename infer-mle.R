#### Libraries and Sources ####

library(ape)
library(diversitree)
library(RPANDA)
library(MLmetrics)
library(ggplot2)
library(bbmle)
source("infer-general-functions.R")


#### end ####


plot_mle_predictions <- function(trees, vec.true.lambda, vec.true.mu, 
                                 lambda_range, mu_range, type = "all", 
                                 save = FALSE){
  
  n_trees <- length(trees)
  n_taxa  <- length(trees[[1]]$tip.label)
  
  vec.pred.rpanda.lambda      <- c()
  vec.pred.rpanda.mu          <- c()
  vec.pred.diversitree.lambda <- c()
  vec.pred.diversitree.mu     <- c()
  vec.pred.ape.lambda         <- c()
  vec.pred.ape.mu             <- c()

  
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

    
  }
  
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
  
  path <- "figures/mle/"
  fname <- paste("mle", type, "ntaxa", n_taxa,
                 "lambda", lambda_range[1], lambda_range[2], 
                 "mu", mu_range[1], mu_range[2],
                 "ntest", n_trees, sep="-")
  fname <- paste(path, fname, sep="")
  
  plot_pred_vs_true(pred.list, true.list, names, "mle", save = save, fname = fname, 
                    subtitle = subtitles)
  
}

