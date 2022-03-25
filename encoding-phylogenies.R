source("infer-general-functions.R")
source("encode-cblv.R")

n_trees       <- 10008    # number of trees to generate
n_taxa        <- c(100,1000) # range size of the generated phylogenies
lambda_range  <- c(0.1, 1.0) # speciation rate
epsilon_range <- c(0.0, 0.9) # turnover rate 
param.range   <- list("lambda"  = lambda_range,
                      "epsilon" = epsilon_range)

out   <- readPhylogeny(n_trees, n_taxa, param.range)
trees <- out$trees
enc   <- generate_encoding(trees, n_taxa)
fname <- getSaveName(n_trees, n_taxa, param.range)$encode
saveRDS(enc, fname)