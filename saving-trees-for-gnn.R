


# Libraries and Sources 
library(svMisc)
source("infer-general-functions.R")

# Simulations parameters 
n_trees  <- 10000 # number of trees 
n_taxa   <- c(100,1000) # taxa range
param.range <- list("lambda" = c(0.1,1.),
                    "epsilon" = c(0.,.9),
                    "q" = c(0.01,0.1)) # range of each parameters 
ss_check <- TRUE # no NAs in the summary statistics?

# Preparation 
out <- load_dataset_trees(n_trees, n_taxa, param.range, ss_check) # load trees 
trees <- out$trees # extract 
list.df.tree <- list() # where datafram will be stored 


#' Create the edges dataframe of a phylo tree 
#'
#' For a graph with N edges, the dataframe has a shape (N, 2). One row 
#' corresponds to one edge, the two columns give the indexes of the two nodes 
#' linked. 
#'
#' @param tree phylo tree 
#'              
#' @return df (edges data frame)
#' 
#' @export
#' @examples
get_edge_df <- function(tree){
  df <- as.data.frame(tree$edge)
  colnames(df) <- c("node1", "node2")
  return(df)
}


#' Create the nodes dataframe of a phylo tree 
#'
#' For a graph with N and considering M attributes per node, the dataframe has a
#' shape (N, M). One row corresponds to one node, each column corresponds to one 
#' attribute. 
#' Examples of nodes attribute: distance to root, number of children, etc.
#'
#' @param tree phylo tree 
#'              
#' @return df (nodes data frame)
#' 
#' @export
#' @examples
get_node_df <- function(tree){
  n_nodes <- 2*tree$Nnode+1
  name.attr <- c("dist", "depth", "ancestor", "descendant")
  n_attr    <- length(name.attr)
  nodes.attr <- vector(mode = "list", length = n_attr)
  names(nodes.attr) <- name.attr
  nodes.attr$dist <- get_all_distances_to_root(tree, as_edge_count = FALSE)
  nodes.attr$depth <- get_all_distances_to_root(tree, as_edge_count = TRUE)
  for (i in 1:n_nodes){
    ancestor <- length(Ancestors(tree, i, "all"))
    descendant <- length(Descendants(tree, i, "all"))
    nodes.attr$ancestor   <- c(nodes.attr$ancestor, ancestor)
    nodes.attr$descendant <- c(nodes.attr$descendant, descendant)
  }
  df <- as.data.frame(nodes.attr)
  return(df)
}


# Computing the node and edge data frame for each tree and save them to the list
for (i in 1:n_trees){
  progress(i, max.value = 10000, progress.bar = TRUE, init = (i == 1))
  tree <- trees[[i]]
  df.edge <- get_edge_df(tree)
  df.node <- get_node_df(tree)
  df.list <- list("node" = df.node, "edge" = df.edge)
  list.df.tree[[i]] <- df.list
}

# Saving 
fname <- "filename.rds" # save name for file 
saveRDS(list.df.tree, fname) # saving file 
car(paste(fname, " saved.\n"))

