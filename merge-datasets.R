
trees.all  <- list()
param.all  <- list()
mle.all    <- list()
ss.all     <- data.frame()
encode.all <- matrix()



for (n in 10000:10009){
  
  print(n)
  
  # File Names 
  fname <- paste("trees-dataset/ntrees", n, "ntaxa-100-1000-lambda-0.1-1-q-0.01-0.1-sscheck-TRUE", sep = "-")
  fname.trees <- paste(fname, "trees.rds", sep = "-")
  fname.param <- paste(fname, "param.rds", sep = "-")
  fname.mle   <- paste(fname, "mle.rds"  , sep = "-")
  fname.ss    <- paste(fname, "ss.rds"  , sep = "-")
  fname.enc   <- paste(fname, "encode.rds"  , sep = "-")

  print("Reading Files...")
  # Reading Files 
  trees      <- readRDS(fname.trees)
  param.true <- readRDS(fname.param)
  param.mle  <- readRDS(fname.mle)
  ss         <- readRDS(fname.ss)
  encode     <- readRDS(fname.enc)
  
  print("Merging Data...")
  # Merging Data
  trees.all <- append(trees.all, trees[1:10000])
  df.param  <- param.true %>% as.data.frame()
  df.mle    <- param.mle  %>% as.data.frame()
  param.all <- rbind(param.all %>% as.data.frame(), df.param[1:10000,]) %>% as.list()
  mle.all   <- rbind(mle.all   %>% as.data.frame(), df.mle[1:10000,])   %>% as.list()
  ss.all    <- rbind(ss.all, ss[1:10000,])
  encode.all <- data.frame(encode.all, encode[,1:10000]) %>% as.matrix()
  #print(dim(encode.all))
}

encode.all <- encode.all[,2:100001]
ss.all <- df_add_tipstate(ss.all, trees.all)

ind <- sample(1:100000, 500)
plot(param.all$q01[ind], mle.all$q01[ind])



fname <- paste("trees-dataset/ntrees", 100000, "ntaxa-100-1000-lambda-0.1-1-q-0.01-0.1-sscheck-TRUE", sep = "-")
fname.trees <- paste(fname, "trees.rds", sep = "-")
fname.param <- paste(fname, "param.rds", sep = "-")
fname.mle   <- paste(fname, "mle.rds"  , sep = "-")
fname.ss    <- paste(fname, "ss.rds"  , sep = "-")
fname.enc   <- paste(fname, "encode.rds"  , sep = "-")
saveRDS(trees.all, fname.trees, compress = TRUE)
saveRDS(param.all, fname.param, compress = TRUE)
saveRDS(mle.all, fname.mle, compress = TRUE)
saveRDS(ss.all, fname.ss, compress = TRUE)
saveRDS(encode.all, fname.enc, compress = TRUE)
