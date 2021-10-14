library(torch)
library(ggplot2)


get_loss_records <- function(fit_model, n_epochs){
  names <- c("epochs", "loss", "type")
  df <- data.frame(matrix(nrow = 0, ncol = length(names)))
  colnames(df) <- names
  for (i in 1:n_epochs){
    train <- fit_model$records$metrics$train[[i]]$loss
    valid <- fit_model$records$metrics$valid[[i]]$loss
    df[nrow(df)+1, ] <- c(i, train, "train")
    df[nrow(df)+1, ] <- c(i, valid, "valid")
  }
  return(df)
}

plot_loss_records <- function(fit_model, n_epochs){
  df_records <- get_loss_records(fitted2, n_epochs)
  
  df_records$epochs <- as.numeric(as.vector(df_records$epochs))
  df_records$loss <- as.numeric(as.vector(df_records$loss))
  ggplot(data=df_records, aes(x=epochs, y=loss, group=type, color=type)) +
    geom_line() + geom_point()+
    scale_y_continuous(trans='log10') +
    theme_minimal()
}