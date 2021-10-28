library(diversitree)
source("infer-general-functions.R")

# Parameteres 
lambda_range  <- c(0.1,1.0)
epsilon_range <- c(0.0,0.9)
n_trees  <- 100000
n_taxa   <- 100
ss_check <- TRUE

# Generate trees 
out <- generate_trees(n_trees, n_taxa, lambda_range, epsilon_range,
                      ss_check = ss_check)


# Extract output data 
trees           <- out$trees
vec.true.lambda <- out$lambda
vec.true.mu     <- out$mu


# Prepare saving - filenames 
save_dataset(trees, vec.true.lambda, vec.true.mu,
             n_trees, n_taxa, lambda_range, epsilon_range, 
             ss_check)

