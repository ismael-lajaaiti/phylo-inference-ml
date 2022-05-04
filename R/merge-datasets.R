source("infer-general-functions.R")

trees.all  <- list()
param.all  <- list()
mle.all    <- list()
ss.all     <- data.frame()
encode.all <- matrix()
df.all     <- list()


for (n in 10000:10009){
  
  print(n)
  
  # File Names 
  fname.list  <- getSaveName(n, n_taxa, param.range)
  fname.trees <- fname.list$trees 
  fname.param <- fname.list$param
  fname.mle   <- fname.list$mle
  fname.ss    <- fname.list$ss
  fname.enc   <- fname.list$encode
  fname.df    <- fname.list$df

  # Reading Files 
  print("Reading Files...")
  trees      <- readRDS(fname.trees)
  param.true <- readRDS(fname.param)
  param.mle  <- readRDS(fname.mle)
  ss         <- readRDS(fname.ss)
  encode     <- readRDS(fname.enc)
  df         <- readRDS(fname.df)
  
  # Merging Data
  print("Merging Data...")
  trees.all  <- append(trees.all, trees[1:10000])
  df.param   <- param.true %>% as.data.frame()
  df.mle     <- param.mle  %>% as.data.frame()
  param.all  <- rbind(param.all %>% as.data.frame(), df.param[1:10000,]) %>% as.list()
  mle.all    <- rbind(mle.all   %>% as.data.frame(), df.mle[1:10000,])   %>% as.list()
  ss.all     <- rbind(ss.all, ss[1:10000,])
  encode.all <- data.frame(encode.all, encode[,1:10000]) %>% as.matrix()
  df.all     <- append(df.all, df[1:10000])
}

encode.all <- encode.all[,2:100001]
ss.all <- df_add_tipstate(ss.all, trees.all)

dim(encode.all)

ind <- sample(1:100000, 500)
plot(param.all$mu[ind], mle.all$mu[ind])



fname.list <- getSaveName(100000, n_taxa, param.range)
fname.trees <- fname.list$trees 
fname.param <- fname.list$param
fname.mle   <- fname.list$mle
fname.ss    <- fname.list$ss
fname.enc   <- fname.list$encode
fname.df    <- fname.list$df
saveRDS(trees.all , fname.trees, compress = TRUE)
saveRDS(param.all , fname.param, compress = TRUE)
saveRDS(mle.all   , fname.mle  , compress = TRUE)
saveRDS(ss.all    , fname.ss   , compress = TRUE)
saveRDS(encode.all, fname.enc  , compress = TRUE)
saveRDS(df.all    , fname.df   , compress = TRUE)

