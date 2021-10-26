# This file contains the general functions needed during the interference 
# processes 


#### Libraries & Sources ####


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
  
  
  if (save){
    fname <- paste(fname, "pdf", sep=".")
    pdf(fname)
    }
  
  n <- length(pred.list)
  n_row <- n %/% 2
  n_col <- 2
  
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





