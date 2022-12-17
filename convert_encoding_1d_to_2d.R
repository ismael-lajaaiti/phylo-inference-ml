source("infer-general-functions.R")

enc_1d <- readRDS("trees-dataset/ntrees-1e+05-ntaxa-100-1000-lambda-0.1-1-q-0.01-0.1-sscheck-TRUE-encode.rds")
concatenate_vector <- c(enc_1d)
enc_2d <- matrix(concatenate_vector, 3*ncol(enc_1d), 1000, byrow=TRUE)
fname <- "trees-dataset/ntrees-1e+05-ntaxa-100-1000-lambda-0.1-1-q-0.01-0.1-sscheck-TRUE-encode2d.rds"
saveRDS(enc_2d, fname)
