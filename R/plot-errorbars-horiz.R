##########################
# SCRIPT TO PLOT FIGURES #
##########################

source("R/infer-phylo-ml.R")

param.crbd      <- readRDS("data/predictions-crbd.rds")
true.param.crbd <- param.crbd$true
pred.param.crbd <- param.crbd$pred
names(pred.param.crbd) <- c("MLE", "CNN-CBLV", "RNN-LTT", "CNN-LTT", "DNN-SS",
                            "GNN-PHY")

param.bisse     <- readRDS("data/predictions-bisse.rds")
true.param.bisse <- param.bisse$true
pred.param.bisse <- param.bisse$pred
names(pred.param.bisse) <- c("MLE", "RNN-LTT", "CNN-LTT", "CNN-CBLV", "DNN-SS",
                            "GNN-PHY")

name.param.crbd  <- names(true.param.crbd)   # name of the predicted parameters
name.method.crbd <- names(pred.param.crbd)   # name of the prediction method
n_param.crbd     <- length(name.param.crbd)  # number of predicted parameters
n_method.crbd    <- length(name.method.crbd) # number of prediction methods

name.param.bisse  <- names(true.param.bisse)   # name of the predicted parameters
name.method.bisse <- names(pred.param.bisse)   # name of the prediction method
n_param.bisse     <- length(name.param.bisse)  # number of predicted parameters
n_method.bisse    <- length(name.method.bisse) # number of prediction methods

range.in.crbd  <- list("lambda"  = c(.2,.9), "mu" = c(0.,.8))
range.in.bisse <- list("lambda0" = c(.2,.9), "q"  = c(0.02,0.09))

error.list.crbd         <- vector(mode = "list", length = n_param.crbd) 
names(error.list.crbd)  <- name.param.crbd
error.list.bisse        <- vector(mode = "list", length = n_param.bisse) 
names(error.list.bisse) <- name.param.bisse

error.list.crbd <- fillErrorList(error.list.crbd, true.param.crbd, 
                                 pred.param.crbd, name.param.crbd,
                                 name.method.crbd, range.in.crbd)

error.list.bisse <- fillErrorList(error.list.bisse, true.param.bisse, 
                                  pred.param.bisse, name.param.bisse,
                                  name.method.bisse, range.in.bisse)


fillErrorList <- function(error.list, true.param, pred.param,
                          name.param, name.method, range.in){
  
  # Set up
  n_param  <- length(name.param)
  n_method <- length(name.method)
  
  for (i in 1:n_param){
    param <- name.param[[i]] # parameter name
    true  <- true.param[[i]] # true values
    range <- range.in[[i]]   # range to study
    for (j in 1:n_method){
      method <- name.method[[j]] # method name
      pred   <- pred.param[[j]][[i]] # predicted values 
      pred_true.filter <- filter_preds_boundary(pred, true, range[1], range[2])
      pred.in <- pred_true.filter$pred.in # predictions in the inner domain
      true.in <- pred_true.filter$true.in # true values in the inner domain 
      theil_coef <- getTheilCoefs(pred.in, true.in) # get the split error 
      error.list[[i]][[method]] <- theil_coef # save split error
    }
  }
  return(error.list)
}


name.error <- c("uniform bias", "consistency bias", "variance", "SSD", "CI")
n_error    <- length(name.method)

df.lambda  <- createEmptyDataFrame(name.error, name.method.crbd)
df.mu      <- createEmptyDataFrame(name.error, name.method.crbd)
df.lambda0 <- createEmptyDataFrame(name.error, name.method.bisse)
df.q       <- createEmptyDataFrame(name.error, name.method.bisse)

createEmptyDataFrame <- function(name.row, name.col){
  n_row <- length(name.row)
  n_col <- length(name.col)
  df    <- data.frame(matrix(ncol = n_col, nrow = n_row))
  colnames(df) <- name.col
  rownames(df) <- name.row
  return(df)
}


fillColumn <- function(df, column, error.list){
  
  df["uniform bias"    , column] <- error.list[[column]]$bias
  df["consistency bias", column] <- error.list[[column]]$slope
  df["variance"        , column] <- error.list[[column]]$var
  df["SSD"             , column] <- error.list[[column]]$SSD
  df["CI"              , column] <- error.list[[column]]$ci
  
  return(df)
  
}

for (method in name.method.crbd){
  df.lambda <- fillColumn(df.lambda, method, error.list.crbd$lambda)
  df.mu     <- fillColumn(df.mu    , method, error.list.crbd$mu)
}

for (method in name.method.bisse){
  df.lambda0 <- fillColumn(df.lambda0, method, error.list.bisse$lambda0)
  df.q       <- fillColumn(df.q      , method, error.list.bisse$q)
}

sortDataFrame <- function(df){
  df <- df %>% t() %>% as.data.frame()
  df <- df %>% arrange(desc(SSD))
  df <- df %>% t() %>% as.data.frame()
  return(df)
}

df.lambda  <- sortDataFrame(df.lambda)
df.mu      <- sortDataFrame(df.mu)
df.lambda0 <- sortDataFrame(df.lambda0)
df.q       <- sortDataFrame(df.q)


legend    <- c("uniform bias", "consistency bias", "variance")
rows <- c("uniform bias", "consistency bias", "variance")

pdf("figures/fig-error-bars-main.pdf", width=6*1.62, height=6)                  # Apply pdf function
par(mfrow = c(2,2),  mar = c(4.5,4,1,1))

plotBar <- function(df, name.param, ylab = "", legend = FALSE, mle_line = FALSE, 
                    xlim = .011){
    ze_barplot <- barplot(as.matrix(df),
                          col = c("black", "darkgray","white"),
                          xlab = "", ylab = "",
                          names.arg = c("", "", "", "", "", ""), xlim = c(0.,xlim),
                          horiz = TRUE, space = .4, angle = 45, border = "black")
    
    grid(nx = NULL, ny = 0, col = "gray", lty = "solid")
    
    ze_barplot <- barplot(as.matrix(df),
                          col = c("#636363", "#D8D8D8","white"),
                          xlab = "", ylab = "",
                          names.arg = c("", "", "", "", "", ""),
                          xlim = c(0., xlim),
                          horiz = TRUE, space = .4, angle = 45, add = T)
    
    
    if (mle_line){abline(v=df["SSD", "MLE"], col="black", lty = "dashed")}
    title(ylab=ylab, line=1, cex.lab=1.7)
    
    xlab <- TeX(paste("$\\", name.param, "$ total error", sep=""))
    title(xlab=xlab, line = 2.5, cex.lab = 1.4)
    
    if (legend){
      legend("topright",                                   
             legend = c("uni. bias", "cons. bias", "variance"),
             fill = c("#636363", "#D8D8D8","white"), 
             inset = c(.0,.0), bg = "white", cex=1.4) 
    }
}

addBarName <- function(df, xpos, inout, space = .4){
  names <- colnames(df)
  n     <- length(names)
  for (i in 1:n){
    y <- 0.9 + (1+space)*(i-1)
    if (inout[i] == "in"){
      text(x = xpos[i], y = y, cex = 1.4, names[i], adj = 1)
    } else if (inout[i] == "out"){
      legend(x = xpos[i], y, names[i], box.col = "white", bg = "white", cex=1.4,
             adj = 0.0, x.intersp=-0.6, y.intersp=-0.2,
             xjust = 0, yjust=0.5, title.adj = 0)
    }
  }
}

addErrorBar <- function(df, space = .4, length = .05){
  names <- colnames(df)
  n     <- length(names)
  for (i in 1:n){
    y        <- 0.9 + (1+space)*(i-1)
    x.center <- df["SSD", i]
    x.width  <- df["CI" , i]
    arrows(x.center - x.width, y, x.center + x.width,
           angle = 90, code = 3, length = length)
  }
}


plotBar(df.lambda[rows,], "lambda", ylab = "CRBD", xlim = .012, legend = TRUE)
xpos.lambda <- as.numeric(df.lambda["SSD",]) - 0.0008
inout.lambda <- c("in", "in", "in", "in", "in", "in")
addBarName(df.lambda, xpos = xpos.lambda, inout = inout.lambda)
addErrorBar(df.lambda)

xpos.mu <- as.numeric(df.mu["SSD",]) + 0.007
xpos.mu[1] <- .123
inout.mu <- c("in", "out", "out", "out", "out", "out")
plotBar(df.mu[rows,], "mu", ylab = "", xlim = .14)
addBarName(df.mu, xpos.mu, inout = inout.mu)
addErrorBar(df.mu)

plotBar(df.lambda0[rows,], "lambda_0", ylab = "BiSSE", xlim = .025)
xpos.lambda0    <- as.numeric(df.lambda0["SSD",]) + 0.0006
xpos.lambda0[1] <- .018
inout.lambda0   <- c("in", "out", "out", "out", "out", "out")
addBarName(df.lambda0, xpos.lambda0, inout = inout.lambda0)
addErrorBar(df.lambda0)

plotBar(df.q[rows,], "q", ylab = "", xlim = .14)
xpos.q    <- as.numeric(df.q["SSD",]) - 0.01
xpos.q[4:6] <- as.numeric(df.q["SSD",4:6]) + 0.003
inout.q   <- c("in", "in", "in", "out", "out", "out")
addBarName(df.q, xpos.q, inout = inout.q)
addErrorBar(df.q)


dev.off()











