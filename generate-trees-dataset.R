source("infer-general-functions.R")


type     <- "bisse"
n_trees  <- 10000 # number of trees to generate
n_taxa   <- c(100,1000) 
ss_check <- TRUE # check that all summary statistics can be computed 
mle      <- TRUE # should mle predictions be computed and saved 


# Constant Rate Birth Death Model 
# -------------------------------
if (type == "crbd"){
  lambda_range  <- c(0.1,1.0)
  epsilon_range <- c(0.0,0.9)
  param.range   <- list("lambda"  = lambda_range,
                        "epsilon" = epsilon_range)
}
# -------------------------------


# Birth Death Model : Exponential Decay Lambda, Constant Mu
# -------------------------------
if (type == "expbd"){
  a_range <- c(.1, .5)
  b_range <- c(.1, .5)
  c_range <- c(.1, .5)
  epsilon_range <- c(0, .9)
  param.range <- list("a"       = a_range,
                      "b"       = b_range,
                      "c"       = c_range,
                      "epsilon" = epsilon_range)
}
# -------------------------------

# BiSSE
# -------------------------------
if (type == "bisse"){
  param.range <- list("lambda" = c(0.1,1.), "epsilon" = c(0.,.9), "q" = c(0.01,0.1))
}
# -------------------------------


# Generate trees 
if(type %in% c("crbd", "expbd")){
  out <- generate_trees(type, n_trees, n_taxa, param.range, ss_check)
}
  
if (type == "bisse"){
  out <- generate_trees_bisse(n_trees, n_taxa, param.range, ss_check)
}
   

# Compute and save tree summary statistics 
if (ss_check){
  df.ss <- generate_ss_dataframe_from_trees(out$trees, out$param)
  save_dataset_summary_statistics(df.ss, n_trees, n_taxa, param.range)
}

# Prepare saving - filenames 
save_dataset_trees(out$trees, out$param, n_trees, n_taxa, param.range, ss_check)

if (mle){
  pred.param <- get_mle_predictions(type, out$trees)
  fname.mle  <- get_mle_preds_save_name(n_trees, n_taxa, param.range, ss_check)
  saveRDS(pred.param, fname.mle)
}


# n_trees
# out <- load_dataset_trees(n_trees, n_taxa, param.range, ss_check)
# true.param <- out$param
# fname.mle  <- get_mle_preds_save_name(n_trees, n_taxa, param.range, ss_check)
# pred.param <- readRDS(fname.mle)
# 
# n_param <- length(true.param)
# par(mfrow = c(1, n_param))
# 
# for (i in 1:n_param){
#   plot(true.param[[i]], pred.param[[i]], ylim = c(0.,0.6))
#   abline(0,1)
# }



