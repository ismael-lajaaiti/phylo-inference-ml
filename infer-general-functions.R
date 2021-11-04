# This file contains the general functions needed during the interference 
# processes 


#### Libraries & Sources ####

library(ape)
library(diversitree)
library(RPANDA)
library(MLmetrics)
library(latex2exp)
source("summary-statistics.R")

#### end ####

#### Operations on Rates ####

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


#### end ####


#### Plotting ####

#' Plot predictions vs. truth
#'
#' Draw a scatter plot represented the predicted values by a model (e.g. DNN) 
#' vs. the true values.
#'
#' @param pred.list list of vectors of predicted paremeters 
#'                  $lambda is the vector of predicted values for lambda
#'                  $mu is the vector of predicted values for mu
#' @param true.list list of vectors of corresponding true values 
#'                  $lambda is the vector of true values for lambda
#'                  $mu is the vector of true values for mu
#' @param save logical, should the plot be saved? (default = FALSE)
#' @param fname if save = TRUE, name of the saved file (default = NA)
#' @param r2_score logical, should the r2 score be printed in the plot? 
#'                 (default = TRUE)
#' @param mse_score logical, should the MSE be printed in the plot?
#'                  (default = TRUE)
#' @param lm_fit logical, should the scatter be fitted by a linear regression?
#'               (default = TRUE)
#' @param subtitles subtitle to print below plots (default = "")
#' @param alone logical, is this function called to draw only one plot? 
#'              (default = TRUE)
#'              
#' @return 
#' 
#' @export
#' @examples
plot_pred_vs_true <- function(pred.list, true.list, names, save = FALSE, fname = NA,
                              r2_score = TRUE, mse_score = TRUE, 
                              lm_fit = TRUE, subtitles = "", alone = TRUE){
  
  n <- length(pred.list)
  n_row <- n %/% 2
  n_col <- 2
  
  lambda.lim <- c(0.05, 1.05) 
  mu.lim <- c(0.,0.9)
  lim <- list(lambda.lim, mu.lim)
  
  if (save){
    fname <- paste(fname, "pdf", sep=".")
    aspect_ratio <- 1.62 * (2 / n_row)
    pdf(fname, width = 10, height = 10/aspect_ratio, pointsize = 15/sqrt(aspect_ratio))
    }
  
  if (alone){par(mfrow=c(n_row,n_col))}
  
  for (i in 1:n){
    
    #pred.name <- paste(names[[1 + (i-1)%%2]], method, sep=".") # parameter name
    pred <- pred.list[[i]] # parameter prediction 
    true <- true.list[[1 + (i-1)%%2]] # parameter true 
    
    # Evaluate the R2 score of predictions vs. truth
    if (r2_score){
      r2 <- R2_Score(pred, true) # compute r2
      r2 <- format(round(r2, 3), nsmall = 3) # format r2 
    }
    
    # Evaluate the MSE of predictions vs. truth 
    if (mse_score){
      mse <- MSE(pred, true) # compute mse 
      mse <- format(round(mse, 3), nsmall = 3) # format mse
    }
    
    rate_name <- names[[1 + (i-1)%%2]]
    main <- create_main_plot_pred_vs_true(rate_name, r2_score, mse_score, r2, mse)
    
    # Scatter plot - Predictions vs. truth
    plot(true, pred, main = main,
         sub = subtitles[1 + (i-1)%/%2], 
         xlim = lim[[1 + (i-1)%%2]], ylim = lim[[1 + (i-1)%%2]])

    # Identity line (for the eye)
    abline(0, 1) # x -> y line (for the eye)
    
    # Linear fit of predictions vs. truth (to see significant trends)
    if (lm_fit){
      fit = lm(pred ~ true)
      fit.scale = lm(pred - true ~ scale(true, scale = FALSE))
      sig = summary(fit)$coefficients[2,4]
      sig.scale = summary(fit.scale)$coefficients[2,4]
      coef <- as.numeric(summary(fit.scale)$coefficients[, 1])
      p_intercept <- as.numeric(summary(fit.scale)$coefficients[1, 4])
      p_slope     <- as.numeric(summary(fit.scale)$coefficients[2, 4])
      sig_intercept <- get_significant_code(p_intercept)
      sig_slope     <- get_significant_code(p_slope)
      intercept <- round(coef[1], 3) 
      slope     <- round(coef[2], 3) 
      text      <- TeX(paste("$\\", slope, "^{", sig_slope, "} \\cdot x + ",
                             intercept, "^{", sig_intercept, sep = ""))
      abline(fit, col="red", lty = ifelse(sig < .05,1,2))
      mtext(text, line = -2, adj = 0)
      #abline(fit.scale, col="blue", lty = ifelse(sig < .05,1,2))
      
    }
  }
  
  if (save){dev.off()}
}


create_main_plot_pred_vs_true <- function(rate_name, r2_score, mse_score,
                                          r2 = NA, mse = NA){
  if (rate_name == "lambda"){main <- '$\\lambda$'}
  else if (rate_name == "mu"){main <- '$\\mu$'}
  if (r2_score){main <- paste(main, '| $r^2 = $', r2)}
  if (mse_score){main <- paste(main, '| $MSE = $', mse)}
  return(TeX(main))
}


get_significant_code <- function(p_value){
  
  s_code <- ""
  if (p_value < 0.001){s_code <- "***"}
  else if (p_value < 0.01){s_code <- "**"}
  else if (p_value < 0.05){s_code <- "*"}
  return(s_code)
  
}



#' Plot the predictions of rates given by MLE
#'
#' Given a list of trees and their true rates, plot the predictions vs. truth
#' where the predictions are obtained by the MLE method.
#'
#' @param trees list of trees 
#' @param vec.true.lambda vector containing the true speciation rates
#' @param vec.true.mu vector containg the true extinction rates 
#' @param lambda_range range in which the speciation rates have been generated
#' @param epsilon_range range in which the turnover rates have been generated 
#' @param type which package to use for the MLE predictions 
#'             either {"ape", "diversitree", "rpanda", "all"}
#'             if "all" use the three packages (ape, diversitree, rpanda) to 
#'             do the predictions
#'             else use only the given package 
#'             (default = "all")
#' @param save logical, should the plot be saved? (default = FALSE)
#' @param fname if save, name of saved file (default = NA)
#' @param alone logical, are the MLE predictions plotted alone or not?
#'              
#' @return 
#' 
#' @export
#' @examples
plot_mle_predictions <- function(trees, vec.true.lambda, vec.true.mu, 
                                 lambda_range, epsilon_range, type = "all", 
                                 save = FALSE, alone = TRUE){
  
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
    lambda_par <- c(.5)  # initial value for lambda
    mu_par     <- c(.4)  # initial value for mu
    
    # RPANDA - Infer lambda & mu 
    if (type == "all" | type == "rpanda"){
      fit.rpanda <- RPANDA::fit_bd(tree, height, f.lambda, f.mu, lambda_par, mu_par, 
                                   cst.lamb = TRUE, cst.mu = TRUE) # fit bd model 
      pred.rpanda.lambda <- fit.rpanda$lamb_par # get lambda prediction 
      pred.rpanda.mu     <- abs(fit.rpanda$mu_par) # get mu prediction 
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
  else{subtitles <- list(paste("MLE:", type), paste("MLE:", type))}
  
  if (type == "rpanda"){
    pred.list <- list(vec.pred.rpanda.lambda, vec.pred.rpanda.mu)
  }
  
  else if (type == "diversitree"){
    pred.list <- list(vec.pred.diversitree.lambda, vec.pred.diversitree.mu)
  }
  
  else if (type == "ape"){
    pred.list <- list(vec.pred.ape.lambda, vec.pred.ape.mu)
  }
  
  else if (type == "all"){
    pred.list <- list(vec.pred.rpanda.lambda, vec.pred.rpanda.mu,
                      vec.pred.diversitree.lambda, vec.pred.diversitree.mu, 
                      vec.pred.ape.lambda, vec.pred.ape.mu)
  }
  
  else{print("Error: type unkown. Type should be either: 'rpanda', 'diversitree', 'ape' or 'all'.")}
  
  
  if (save){
    fname <- get_plot_save_name("mle", n_trees, n_taxa, lambda_range, epsilon_range,
                                n_test)
  }
  
  plot_pred_vs_true(pred.list, true.list, names, save = save, fname = fname, 
                    subtitle = subtitles, alone = alone)
  
}


#' Plot together predictions from Neural Networks and MLE 
#' 
#'
#' @param pred.nn.list list of predicted parameters by the Neural Network
#' @param true.list true values of the parameters 
#' @param trees list of trees on which the predictions of rates have been made
#' @param n_trees number of trees in trees 
#' @param n_taxa size of the trees 
#' @param lambda_range range in which the speciation rates have been generated
#' @param epsilon_range range in which the turnover rates have been generated 
#' @param mle_package which MLE package to use (either "ape", "diversitree", "rpanda")
#'                    (default = "ape")
#' @param nn_type type of Neural Network (e.g. CNN, DNN, RNN) (default = "NN")
#' @param save logical, should the plot be saved? (default = FALSE)
#' @param fname if save, file name
#' 
#'              
#' @return fname, file name of the prediction pots 
#' 
#' @export
#' @examples
plot_together_nn_mle_predictions <- function(pred.nn.list, true.list, trees, nn_type,
                                             n_trees, n_taxa, lambda_range, epsilon_range,
                                             n_layer, n_hidden, n_train, ker_size = NA,
                                             mle_package = "ape", 
                                             save = FALSE){
  
  if (save){
    n_test <- length(trees)
    fname <- get_plot_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range,
                                n_test, n_layer, n_hidden, n_train, ker_size)
    aspect_ratio <- 1.62
    pdf(fname, width = 10, height = 10/aspect_ratio, pointsize = 15/sqrt(aspect_ratio))
  }
  
  name.list <- list("lambda", "mu")
  
  par(mfrow=c(2,2))
  plot_pred_vs_true(pred.nn.list, true.list, name.list,
                    subtitles = c(nn_type, nn_type), save = FALSE,
                    alone = FALSE)
  vec.true.lambda <- true.list[[1]]
  vec.true.mu     <- true.list[[2]]
  plot_mle_predictions(trees, vec.true.lambda, vec.true.mu, 
                       lambda_range, epsilon_range, 
                       type = mle_package, save = FALSE, alone = FALSE)
  
  if (save){dev.off()}
  
}


#### end ####


#### Generating, Converting, Saving and Loading ####

#' Generate phylogenetic trees 
#'
#' Create a given number of trees of a choosen size and rates. This trees are
#' stored in a list.
#'
#' @param n_trees number of trees to generate 
#' @param n_taxa size of the trees 
#' @param lambda_range vector containg the min and max of the speciation rates 
#'                     of the trees generated. For each trees the speciation 
#'                     rate will randomly drawn in this interval
#' @param lambda_range vector containg the min and max of the turnover rates 
#'                     of the trees generated. For each trees the turnover 
#'                     rate will randomly drawn in this interval
#' @param ss_check logical, should we check that the summary statistics 
#'                 corresponding to the generated trees doesn't contain any NA
#'                 values (default = TRUE)
#'              
#' @return out list of outputs 
#'         $trees -> trees, the list of the generated trees 
#'         $lambda -> vec.lambda, vector of speciation rates of the generated trees
#'         $mu -> vec.mu, vector of extinction rates of generated trees
#' 
#' @export
#' @examples
generate_trees <- function(n_trees, n_taxa, lambda_range, epsilon_range,
                           ss_check = TRUE){
  
  trees <- list() # initialize tree list where trees will be stored 
  vec.true.lambda <- c() # vector where true lambda values will be stored
  vec.true.mu     <- c() # vector where true mu values will be stored 
  
  print(n_taxa)
  print(length(n_taxa))

  cat("Generation of trees...\n")
  
  while (length(trees) < n_trees){
    # Generate the phylogenetic tree 
    true.lambda  <- runif(1, lambda_range[1] , lambda_range[2]) # generate random lambda
    epsilon      <- runif(1, epsilon_range[1], epsilon_range[2]) # generate random epsilon
    true.mu      <- epsilon * true.lambda # compute corresponding mu
    if (length(n_taxa) == 1){
      tree <- trees(c(true.lambda, true.mu), "bd", max.taxa=n_taxa)[[1]] # create the tree 
    }
    else if (length(n_taxa) == 2) {
      true.n_taxa <- sample(n_taxa[1]:n_taxa[2], 1)
      tree <- trees(c(true.lambda, true.mu), "bd",
                    max.taxa=true.n_taxa)[[1]] # create the tree
    }
    else (cat("Error n_taxa should either be an integer, or a vector of size 2"))
    
    
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


#' Generate a data.frame containing Lineage Through Time of phylo trees
#'
#' Given a list of phylo trees, returns their LTT in a data.frame
#'
#' @param trees list of phylo trees
#' @param n_taxa size of the trees 
#' @param vec.lambda vector of the corresponding speciation rates of the
#'                   generated trees 
#' @param vec.mu vector of the corresponding extinction rates of the
#'                   generated trees                   
#'              
#' @return out list of outputs
#'         $ltt -> df.ltt data.frame containing the LTT of the trees 
#'                i'th column contains the LTT of the i'th tree
#'         $rates -> df.rates data.frame containing the speciation and extinction rates  
#' 
#' @export
#' @examples
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


#' Convert a LTT data.frame into a torch dataset 
#'
#' Given LTT data.frame and a data.frame containing rates (targets to infer)
#' create the corresponding torch dataset 
#'
#' @param df.ltt data.frame containing LTT (created by generate_ltt_dataframe)
#' @param df.rates data.frame containing rates (created by generated_ltt_dataframe)
#' @param vec.lambda vector of the corresponding speciation rates of the
#'              
#' @return torch dataset 
#'         $x -> LTT contained in df.ltt
#'         $y -> corresponding rates contained in df.rates (targets)
#' 
#' @export
#' @examples
convert_ltt_dataframe_to_dataset <- function(df.ltt, df.rates){
  
  ds.ltt <- torch::dataset(
    name <- "ltt_dataset", 
    initialize = function(df.ltt, df.rates){
      
      array.ltt <- df.ltt %>% 
        as.matrix() %>% 
        array(dim = c(nrow(df.ltt), ncol(df.ltt), 1)) # convert df to array of good dimension
      
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
      list(x = self$x[,i ,], y = self$y[i, ])
    }, 
    
    .length = function() {
      self$y$size()[[1]]
    }
  )
  
  return(ds.ltt)
}



convert_ltt_dataframe_to_dataset_cnn <- function(df.ltt, df.rates){
  
  ds.ltt <- torch::dataset(
    name <- "ltt_dataset", 
    initialize = function(df.ltt, df.rates){
      
      array.ltt <- df.ltt %>% 
        as.matrix() %>% 
        torch_tensor() # convert df to array of good dimension
      
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
      list(x = self$x[ , i]$unsqueeze(1), y = self$y[i, ])
    }, 
    
    .length = function() {
      self$y$size()[[1]]
    }
  )
  
  return(ds.ltt)
}

get_backbone_save_name <- function(n_trees, n_taxa, lambda_range, epsilon_range){
  
  
  if (length(n_taxa) == 1){
    fname <- paste("ntrees", n_trees, "ntaxa", n_taxa,
                   "lambda" , lambda_range[1] , lambda_range[2], 
                   "epsilon", epsilon_range[1], epsilon_range[2], sep="-")
  }
  
  else if (length(n_taxa) == 2){
    fname <- paste("ntrees", n_trees, "ntaxa", n_taxa[1], n_taxa[2],
                   "lambda" , lambda_range[1] , lambda_range[2], 
                   "epsilon", epsilon_range[1], epsilon_range[2], sep="-")
  }
  
  else {cat("Error length(n_taxa) should be either equal to 1 or 2. ")}
  
  return(fname)
  
}


get_model_save_name <- function(nn_type, n_trees, n_taxa, lambda_range, epsilon_range,
                                 n_layer, n_hidden, n_train, ker_size = NA){
  
  fname <- get_backbone_save_name(n_trees, n_taxa, lambda_range, epsilon_range)
  fname.model <- paste(fname, "nlayer", n_layer, "nhidden", n_hidden, 
                       "ntrain", n_train, sep = "-")
  if (!is.na(ker_size)){fname.model <- paste(fname.model, "kersize", ker_size, sep = "-")}
  dir <- paste("neural-networks-models", nn_type, "", sep = "/")
  fname.model <- paste(dir, fname.model, sep = "")
  return(fname.model)
}


get_dataset_save_name <- function(n_trees, n_taxa, lambda_range, epsilon_range,
                                  ss_check){
  
  dir <- "trees-dataset/"
  
  fname <- get_backbone_save_name(n_trees, n_taxa, lambda_range, epsilon_range)
  fname <- paste(fname, "sscheck", ss_check, sep = "-")
  
  fname.trees  <- paste(dir, fname, "-trees.rds", sep="")
  fname.lambda <- paste(dir, fname, "-lambda.rds", sep="")
  fname.mu     <- paste(dir, fname, "-mu.rds", sep="")
  fname.ss     <- paste(dir, fname, "-ss.rds", sep="")
  
  fnames <- list("trees"  = fname.trees, 
                 "lambda" = fname.lambda, 
                 "mu"     = fname.mu,
                 "ss"     = fname.ss)
  
  return(fnames)
  
}


#' Save generated trees and their rates 
#'
#' Given a list of trees and their rates, those are saved using saveRDS function
#'
#' @param trees list of trees 
#' @param vec.lambda vector containing the speciation rates
#' @param vec.mu vector containg the extinction rates 
#' @param n_trees number of trees in trees 
#' @param n_taxa size of the trees 
#' @param lambda_range range in which the speciation rates have been generated
#' @param epsilon_range range in which the turnover rates have been generated 
#' @param ss_check logical, have we check that the trees generated doesn't 
#'                 give NA values for the summary statistics
#'              
#' @return 
#' 
#' @export
#' @examples
save_dataset_trees <- function(trees, vec.lambda, vec.mu, n_trees, n_taxa,
                         lambda_range, epsilon_range, ss_check){
  
  # Getting file names to save 
  fnames <- get_dataset_save_names(n_trees, n_taxa, lambda_range, epsilon_range,
                                   ss_check)
  fname.trees  <- fnames$trees
  fname.lambda <- fnames$lambda
  fname.mu     <- fnames$mu
  
  # Saving data 
  cat("Saving data...\n")
  saveRDS(trees, fname.trees)
  cat(paste(fname.trees, " saved.\n", sep=""))
  saveRDS(vec.lambda, fname.lambda)
  cat(paste(fname.lambda, " saved.\n", sep=""))
  saveRDS(vec.mu, fname.mu)
  cat(paste(fname.mu, " saved.\n", sep=""))
  cat("Saving data... Done.\n")
  
}


#' Load trees list and rates vector
#'
#' Load a saved list of trees and corresponding rates, 
#' given the information needed to recover their names
#'
#' @param n_trees number of trees in trees 
#' @param n_taxa size of the trees 
#' @param lambda_range range in which the speciation rates have been generated
#' @param epsilon_range range in which the turnover rates have been generated 
#' @param ss_check ss_check logical, have we check that the trees generated doesn't 
#'                 give NA values for the summary statistics
#'              
#' @return out list of outputs 
#'         $trees -> trees, the list of the generated trees 
#'         $lambda -> vec.lambda, vector of speciation rates of the generated trees
#'         $mu -> vec.mu, vector of extinction rates of generated trees
#' 
#' @export
#' @examples
load_dataset_trees <- function(n_trees, n_taxa, lambda_range, mu_range, ss_check){
  
  # Getting file names to load 
  fnames <- get_dataset_save_name(n_trees, n_taxa, lambda_range, mu_range,
                                   ss_check)
  fname.trees  <- fnames$trees
  fname.lambda <- fnames$lambda
  fname.mu     <- fnames$mu
  
  # Loading data 
  cat("Loading data...\n")
  trees <- readRDS(fname.trees)
  cat(paste(fname.trees, " loaded.\n", sep=""))
  vec.lambda <- readRDS(fname.lambda)
  cat(paste(fname.lambda, " loaded.\n", sep=""))
  vec.mu <- readRDS(fname.mu)
  cat(paste(fname.mu, " loaded.\n", sep=""))
  cat("Loading data... Done.\n")
  
  out <- list("trees" = trees, "lambda" = vec.lambda, "mu" = vec.mu)
  
  return(out)
  
}


#' Save computed trees summary statistics data.frame
#'
#' 
#'
#' @param df.ss data.frame containing the summary statics and rates of several trees
#' @param n_trees number of trees in trees 
#' @param n_taxa size of the trees 
#' @param lambda_range range in which the speciation rates have been generated
#' @param epsilon_range range in which the turnover rates have been generated 
#'              
#' @return 
#' 
#' @export
#' @examples
save_dataset_summary_statistics <- function(df.ss, n_trees, n_taxa,
                                            lambda_range, epsilon_range){
  
  fnames <- get_dataset_save_name(n_trees, n_taxa, lambda_range, epsilon_range,
                                  ss_check = TRUE)
  fname.ss <- fnames$ss
  
  cat("Saving summary statistics data...\n")
  saveRDS(df.ss, fname.ss)
  cat(paste(fname.ss, " saved.\n", sep=""))
  cat("Saving summary statistics data... Done.\n")
  
}


#' Load a saved summary statistics data.frame
#'
#' Load a summary statistics data.frame, given the information needed to recover
#' its name 
#'
#' @param n_trees number of trees in trees 
#' @param n_taxa size of the trees 
#' @param lambda_range range in which the speciation rates have been generated
#' @param epsilon_range range in which the turnover rates have been generated 
#'              
#' @return df.ss summary statistics data.frame
#' 
#' @export
#' @examples
load_dataset_summary_statistics <- function(n_trees, n_taxa,
                                            lambda_range, epsilon_range){
  
  fnames <- get_dataset_save_name(n_trees, n_taxa, lambda_range, epsilon_range,
                                  ss_check = TRUE)
  fname.ss <- fnames$ss
  
  cat("Loading summary statistics data...\n")
  df.ss <- readRDS(fname.ss)
  cat("Loading summary statistics data... Done.\n")
  
  return(df.ss)
  
}


save_predictions <- function(pred.list, true.list, nn_type, n_trees, n_taxa,
                             lambda_range, epsilon_range, n_test, n_layer, 
                             n_hidden, n_train, ker_size = NA){
  
  cat("Saving predictions...\n")
  save.list <- list("pred" = pred.list, "true" = true.list)
  fname.preds <- get_preds_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range,
                                          n_test, n_layer, n_hidden, n_train, ker_size)
  saveRDS(save.list, fname.preds)
  cat("\nSaving predictions... Done.\n")
  
}


load_predictions <- function(nn_type, n_trees, n_taxa,
                             lambda_range, epsilon_range, n_test, n_layer, 
                             n_hidden, n_train, ker_size = NA){
  
  fname.preds <- get_preds_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range,
                                     n_test, n_layer, n_hidden, n_train, ker_size)
  
  save.list <- readRDS(fname.preds)
  
  return(save.list)
  
}


#' Create file name for the prediction plots 
#'
#'
#' @param n_trees number of trees in trees 
#' @param n_taxa size of the trees 
#' @param lambda_range range in which the speciation rates have been generated
#' @param epsilon_range range in which the turnover rates have been generated 
#' @param n_test number of tested predictions
#' @param dir directory where the file will be saved 
#' @param type if the model is a neural network, type should be "nn" such that,
#'             relevant informations will be added to the file name
#' @param n_layer if type="nn", number of layer of the neural network (default = NA)
#' @param n_hidden if type="nn", number of neurons in hidden layeres (default = NA)
#' @param n_train if type="nn", size of the training set (default = NA)
#' 
#'              
#' @return fname, file name of the prediction pots 
#' 
#' @export
#' @examples
get_plot_save_name <- function(model_type, n_trees, n_taxa, lambda_range, epsilon_range,
                               n_test, n_layer = NA, n_hidden = NA, n_train = NA, ker_size = NA){
  
  fname <- get_backbone_save_name(n_trees, n_taxa, lambda_range, epsilon_range)

  if (model_type != "mle") {
    fname.plot <- paste(fname, "nlayer", n_layer, "nhidden", n_hidden, 
                       "ntrain", n_train, sep = "-")
    if (!is.na(ker_size)){fname.plot <- paste(fname.plot, "kersize", ker_size, sep = "-")}
  }

  dir.plot <- paste("figures", model_type, "", sep = "/")
  fname.plot <- paste(dir.plot, fname.plot, "-ntest-", n_test, sep = "")
  fname.plot <- paste(fname.plot, "pdf", sep = ".")
  
  return(fname.plot)
  
}


get_preds_save_name <- function(nn_type, n_trees, n_taxa, lambda_range, epsilon_range,
                                 n_test, n_layer, n_hidden, n_train, ker_size = NA){
  
  fname.plot <- get_plot_save_name(nn_type, n_trees, n_taxa, lambda_range, epsilon_range,
                                   n_test, n_layer, n_hidden, n_train, ker_size)
  
  start <- nchar("figures") + nchar(nn_type) + 3
  fname.preds <- substring(fname.plot, start, nchar(fname.plot) - 4)
  dir.preds   <- paste("neural-networks-predictions", nn_type, "", sep = "/")
  fname.preds <- paste(dir.preds, fname.preds, ".rds", sep = "")
  
  return(fname.preds)
}


#### end #### 







