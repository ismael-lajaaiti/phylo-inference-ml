#### File description ####
# Generate phylogenies either with the Constant Rate Birth Death model (CRBD)
# or with the Binary State Speciation and Extinction model (BiSSE)
#### end ####



#### Import libraries #### 

source("infer-general-functions.R")

#### end ####


model       <- "crbd"      # type of the model, either: "crbd" or "bisse"
n_trees     <- 10009      # number of trees to generate
n_taxa      <- c(100,1000) # range size of the generated phylogenies
compute_mle <- TRUE        # should mle predictions be computed and saved 


#### Defining parameter space of the model ####

# For the CRBD model
lambda_range       <- c(0.1, 1.0) # speciation rate
epsilon_range      <- c(0.0, 0.9) # turnover rate 
param.range.crbd   <- list("lambda"  = lambda_range,
                           "epsilon" = epsilon_range)

# For the BiSSE model 
lambda_range      <- c(0.1, 1.)  # speciation rate
q_range           <- c(0.01,0.1) # transition rate 
param.range.bisse <- list("lambda" = lambda_range, 
                          "q"      = q_range)

# Select the parameter space of the choosen diversification model
param.range.list <- list("crbd"  = param.range.crbd,
                         "bisse" = param.range.bisse)
param.range <- param.range.list[[model]] 

#### end ####


#### Generate phylogenies - Compute Sum. Stat. & MLE ####

# Generating and saving phylogenies
out <- generatePhylo(model, n_trees, n_taxa, param.range)
savePhylogeny(out$trees, out$param, n_trees, n_taxa, param.range) # save


# Computing summary statistics 
df.ss <- generateSumStatFromPhylo(out$trees, out$param)      # compute
saveSummaryStatistics(df.ss, n_trees, n_taxa, param.range) # save


# Computing Maximum Likelihood Rate Estimations
if (compute_mle){
  pred.param <- getPredsMLE(model, out$trees)
  savePredsMLE(pred.param, n_trees, n_taxa, param.range)
}

#### end ####
