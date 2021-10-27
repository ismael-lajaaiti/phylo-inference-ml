library(torch)
library(ggplot2)


get_loss_records <- function(net.fit){
  names <- c("epochs", "loss", "type")
  df <- data.frame(matrix(nrow = 0, ncol = length(names)))
  colnames(df) <- names
  n <- length(net.fit$records$metrics$train)  # number of epochs 
  for (i in 1:n){
    train <- net.fit$records$metrics$train[[i]]$loss
    valid <- net.fit$records$metrics$valid[[i]]$loss
    df[nrow(df)+1, ] <- c(i, train, "train")
    df[nrow(df)+1, ] <- c(i, valid, "valid")
  }
  return(df)
}


plot_loss_records <- function(net.fit){
  df_records <- get_loss_records(net.fit)
  
  df_records$epochs <- as.numeric(as.vector(df_records$epochs))
  df_records$loss <- as.numeric(as.vector(df_records$loss))
  ggplot(data=df_records, aes(x=epochs, y=loss, group=type, color=type)) +
    geom_line() + geom_point()+
    scale_y_continuous(trans='log10') +
    theme_minimal()
}






