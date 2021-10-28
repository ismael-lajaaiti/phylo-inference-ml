library(diversitree)
source("infer-general-functions.R")

# Parameteres 
lambda_range  <- c(0.1,1.0)
epsilon_range <- c(0.0,0.9)
n_trees  <- 10000
n_taxa   <- 300
ss_check <- TRUE

# Generate trees 
out <- generate_trees(n_trees, n_taxa, lambda_range, epsilon_range,
                      ss_check = ss_check)


# Extract output data 
trees           <- out$trees
vec.true.lambda <- out$lambda
vec.true.mu     <- out$mu


# Prepare saving - filenames 
fnames <- get_dataset_save_names(n_trees, n_taxa, lambda_range, mu_range, ss_check)
fname.trees  <- fnames$trees
fname.lambda <- fnames$lambda
fname.mu     <- fnames$mu

# Saving data 
cat("Saving data...\n")
saveRDS(trees, fname.trees)
cat(paste(fname.trees, " saved.\n", sep=""))
saveRDS(vec.true.lambda, fname.lambda)
cat(paste(fname.lambda, " saved.\n", sep=""))
saveRDS(vec.true.mu, fname.mu)
cat(paste(fname.mu, " saved.\n", sep=""))
cat("\nSaving data... Done.\n")

