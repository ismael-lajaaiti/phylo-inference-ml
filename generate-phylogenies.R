# Simulate phylogenies either with the Constant Rate Birth Death model (CRBD)
# or with the Binary State Speciation and Extinction model (BiSSE)

#### Dependencies ####

source("R/libraries.R")
source("infer-general-functions.R")

#### end ####

model <- "bisse" # type of the model, either: "crbd" or "bisse"
n_trees <- 100000 # number of trees to generate
n_taxa <- c(100, 1000) # range size of the generated phylogenies
compute_mle <- TRUE # should mle predictions be computed and saved

#### Defining parameter space of the model ####

# For the CRBD model
lambda_range <- c(0.1, 1.0) # speciation rate
epsilon_range <- c(0.0, 0.9) # turnover rate
param.range.crbd <- list(
    "lambda" = lambda_range,
    "epsilon" = epsilon_range
)

# For the BiSSE model
lambda_range <- c(0.1, 1.) # speciation rate
q_range <- c(0.01, 0.1) # transition rate
param.range.bisse <- list(
    "lambda" = lambda_range,
    "q" = q_range
)

# Select the parameter space of the choosen diversification model
param.range.list <- list(
    "crbd" = param.range.crbd,
    "bisse" = param.range.bisse
)
param.range <- param.range.list[[model]]

#### end ####


#### Generate phylogenies - Compute Sum. Stat. & MLE ####

# Generating and saving phylogenies

n_core <- 10
n_rep <- 10
n_trees_per_rep <- 50000

r <- mclapply(1:n_rep, function(i) {
    out <- generatePhylo(model, n_trees_per_rep, n_taxa, param.range)
    print("Phylogenies generated.")
    mle.param <- getPredsMLE(model, out$trees)
    print("MLE predictions computed.")
    list(out = out, mle = mle.param)
}, mc.cores = n_core)

for (i in 1:n_rep) {
    out <- r[[i]]$out
    mle <- r[[i]]$mle
    ss <- out$ss
    trees <- out$trees
    params <- out$param
    trees_params <- list(trees = trees, params = params)
    fname_prefix <- paste("trees-dataset/bisse-n", n_trees_per_rep, sep = "")
    i_pad <- str_pad(i, 2, pad = "0") # ex: 2 -> 02 | 13 -> 13
    saveRDS(trees_params, paste(fname_prefix, "-trees-and-params", i_pad, ".rds", sep = ""))
    saveRDS(ss, paste(fname_prefix, "-sumstat", i_pad, ".rds", sep = ""))
    saveRDS(mle, paste(fname_prefix, "-predmle", i_pad, ".rds", sep = ""))
}



#
# savePhylogeny(out$trees, out$param, n_trees, n_taxa, param.range) # save
#
#
# # Computing summary statistics
# df.ss <- generateSumStatFromPhylo(out$trees, out$param)      # compute
# saveSummaryStatistics(df.ss, n_trees, n_taxa, param.range) # save
#
#
# # Computing Maximum Likelihood Rate Estimations
# if (compute_mle){
#   pred.param <- getPredsMLE(model, out$trees)
#   savePredsMLE(pred.param, n_trees, n_taxa, param.range)
# }
#
# #### end ####
