source("~/phylo-inference-ml/encode-cblv.R")

n <- 10009 # number of trees
n_taxa <- c(100,1000) # taxa range 


# Loading phylo trees to encode
fname <- paste("trees-dataset/ntrees", n, "ntaxa-100-1000-lambda-0.1-1-q-0.01-0.1-sscheck-TRUE", sep = "-")
fname.trees <- paste(fname, "trees.rds", sep = "-")
trees <- readRDS(fname.trees)

# Generate and save the encoding 
encode.matrix <- generate_encoding_bisse(trees, n_taxa) 
fname.encode  <- paste(fname, "encode.rds", sep = "-")
saveRDS(encode.matrix, fname.encode)





