##########################
# SCRIPT TO PLOT FIGURES #
##########################





param      <- readRDS("params-testset-bisse-encode.rds")
true.param <- param$true
pred.param <- param$pred
names(pred.param) <- c("MLE", "CNN-CBLV", "CNN-CBLV-disor", "CNN-CBLV-notips" ,"CNN-LTT")

name.param  <- names(true.param)   # name of the predicted parameters
name.method <- names(pred.param)   # name of the prediction method
n_param     <- length(name.param)  # number of predicted parameters
n_method    <- length(name.method) # number of prediction methods

#range.in.crbd  <- list("lambda"  = c(.2,.9), "mu" = c(0.,.8))
range.in <- list("lambda0" = c(.2,.9), "q"  = c(0.02,0.09))

error.list         <- vector(mode = "list", length = n_param) 
names(error.list)  <- name.param

error.list <- fillErrorList(error.list, true.param, 
                            pred.param, name.param,
                            name.method, range.in)


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
n_error    <- length(name.error)

df.lambda0 <- createEmptyDataFrame(name.error, name.method)
df.q       <- createEmptyDataFrame(name.error, name.method)

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

for (method in name.method){
  df.lambda0 <- fillColumn(df.lambda0, method, error.list$lambda0)
  df.q       <- fillColumn(df.q      , method, error.list$q)
}

sortDataFrame <- function(df){
  df <- df %>% t() %>% as.data.frame()
  df <- df %>% arrange(desc(SSD))
  df <- df %>% t() %>% as.data.frame()
  return(df)
}

df.lambda0 <- sortDataFrame(df.lambda0)
df.q       <- sortDataFrame(df.q)


legend    <- c("uniform bias", "consistency bias", "variance")
rows <- c("uniform bias", "consistency bias", "variance")

pdf("testtt.pdf", width=6*1.62, height=3)                  # Apply pdf function
par(mfrow = c(1,2),  mar = c(4.5,4,1,1))

plotBar <- function(df, name.param, ylab = "", legend = FALSE, mle_line = FALSE, 
                    xlim = .011){
  ze_barplot <- barplot(as.matrix(df),
                        col = c("black", "darkgray","white"),
                        xlab = "", ylab = "",
                        names.arg = c("", "", "", "", ""), xlim = c(0.,xlim),
                        horiz = TRUE, space = .4, angle = 45, border = "black")
  
  grid(nx = NULL, ny = 0, col = "gray", lty = "solid")
  
  ze_barplot <- barplot(as.matrix(df),
                        col = c("#636363", "#D8D8D8","white"),
                        xlab = "", ylab = "",
                        names.arg = c("", "", "", "", ""),
                        xlim = c(0., xlim),
                        horiz = TRUE, space = .4, angle = 45, add = T)
  
  
  if (mle_line){abline(v=df["SSD", "MLE"], col="black", lty = "dashed")}
  title(ylab=ylab, line=1, cex.lab=1.4)
  
  xlab <- TeX(paste("$\\", name.param, "$ total error", sep=""))
  title(xlab=xlab, line = 2.5, cex.lab = 1.)
  
  if (legend){
    legend("topright",                                   
           legend = c("uni. bias", "cons. bias", "variance"),
           fill = c("#636363", "#D8D8D8","white"), 
           inset = c(.0,.0), bg = "white", cex=1.) 
  }
}

addBarName <- function(df, xpos, inout, space = .4){
  names <- colnames(df)
  n     <- length(names)
  for (i in 1:n){
    y <- 0.9 + (1+space)*(i-1)
    name <- names[i]
    if (names[i] == "CNN-CBLV-disor"){
      name <- TeX("CNN-$CBLV^{disotips}$")
      y <- y +.05*as.integer(inout[i] == "in")
      }
    if (
      names[i] == "CNN-CBLV-notips"){name <- TeX("CNN-$CBLV^{notips}$")
      y <- y +.05*as.integer(inout[i] == "in")
      }
    if (inout[i] == "in"){
      text(x = xpos[i], y = y, cex = 1., name, adj = 1.)
    } else if (inout[i] == "out"){
      legend(x = xpos[i], y, name, box.col = "white", bg = "white", cex=1.,
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


plotBar(df.lambda0[rows,], "lambda_0", ylab = "BiSSE", xlim = .012, legend = TRUE)
xpos.lambda0 <- as.numeric(df.lambda0["SSD",]) - 0.0008
xpos.lambda0[3:5] <- as.numeric(df.lambda0["SSD",3:5]) + 0.0005
inout.lambda0 <- c("in", "in", "out", "out", "out")
addBarName(df.lambda0, xpos = xpos.lambda0, inout = inout.lambda0)
addErrorBar(df.lambda0)

plotBar(df.q[rows,], "q", ylab = "", xlim = .14)
xpos.q    <- as.numeric(df.q["SSD",]) - 0.01
xpos.q[3:5] <- as.numeric(df.q["SSD",3:5]) + 0.003
inout.q   <- c("in", "in", "out", "out", "out")
addBarName(df.q, xpos.q, inout = inout.q)
addErrorBar(df.q)


dev.off()