source("infer-general-functions.R")
source("encode-cblv.R")
library(torch)
library(parallel)
library(stringr)


#### Defining parameter space of the model ####

n_core <- 20
n_rep <- 30
n_trees_per_rep <- 20000
n_taxa <- c(100, 1000)

r <- mclapply(1:n_rep, function(i){
  fname_prefix <- paste("data-from-7/bisse-n", n_trees_per_rep, sep="")
  i_pad <- str_pad(i+20, 2, pad="0") # ex: 2 -> 02 | 13 -> 13 
  print(i_pad)
  out <- readRDS(paste(fname_prefix, "-trees-and-params", i_pad,".rds",sep=""))
  trees <- out$trees
  enc <- generate_encoding_bisse(trees, n_taxa)
  print("Phylogenies encoded.")
  list(enc=enc)
}, mc.cores=n_core)

for (i in 1:n_rep){
  enc <- r[[i]]$enc
  fname_prefix <- paste("data-from-7/bisse-n", n_trees_per_rep, sep="")
  i_pad <- str_pad(i+20, 2, pad="0") # ex: 2 -> 02 | 13 -> 13 
  saveRDS(enc, paste(fname_prefix, "-cblv", i_pad,".rds",sep=""))
}