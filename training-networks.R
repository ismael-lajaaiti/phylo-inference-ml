#### Training loop ####


trainingNetwork <- function(n_epochs, patience, training_set, validation_set,
                            model_type){

  if (model_type == "dnn_ss"){trainingDNNss(n_epochs, patience,)}

  else if (model_type == "cnn_ltt"){trainingCNNltt(n_epochs, patience)}

  else if (model_type == "rnn_ltt"){trainingRNNltt(n_epochs, patience)}

  else if (model_type == "cnn_enc"){trainingCNNenc(n_epochs, patience)}

  else{ print("Error: Model Type unkown."}

}

trainingDNNss <- function(n_epochs, patience){



}



#### end ####
