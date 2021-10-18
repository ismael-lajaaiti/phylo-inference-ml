get_r_epsilon <- function(vec.lambda_mu){
  r <- lambda - mu
  epsilon <- mu / lambda
  vec.r_epsilon <- c(r, epsilon)
  return(vec.r_epsilon)
}

get_lambda_mu <- function(vec.r_epsilon){
  lambda <- r / (1 - epsilon)
  mu <- epsilon * lambda 
  vec.lambda_mu <- c(lambda, mu)
  return(vec.lambda_mu)
}