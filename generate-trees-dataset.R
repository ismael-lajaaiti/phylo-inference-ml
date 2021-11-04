library(diversitree)
source("infer-general-functions.R")

# Parameteres 
lambda_range  <- c(0.1,1.0)
epsilon_range <- c(0.0,0.9)
n_trees  <- 10000
n_taxa   <- c(100,400)
ss_check <- TRUE

# Generate trees 
out <- generate_trees(n_trees, n_taxa, lambda_range, epsilon_range,
                      ss_check = ss_check)

# Extract output data 
trees           <- out$trees
vec.true.lambda <- out$lambda
vec.true.mu     <- out$mu

if (ss_check){
  df.ss <- generate_ss_dataframe_from_trees(trees, vec.true.lambda, vec.true.mu)
  save_dataset_summary_statistics(df.ss, n_trees, n_taxa, lambda_range, 
                                  epsilon_range)
}

# Prepare saving - filenames 
save_dataset_trees(trees, vec.true.lambda, vec.true.mu,
                   n_trees, n_taxa, lambda_range, epsilon_range, 
                   ss_check)

