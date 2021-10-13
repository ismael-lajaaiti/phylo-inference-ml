# This file contains the function to encode a phylogenetic tree with the
# "Compact Bijective Ladderized Vector" methods
# See Voznica 2021 - bioRxiv - DOI:10.1101/2021.03.11.435006 


#### Libraries ####


library(ape)
library(diversitree)
library(castor)
library(phangorn)


#### end ####


#### Encoding Tree ####

#' Check a node of a tree is a tip
#'
#'
#' @param node integer, index of the node
#' @param tree phylo tree
#'
#' @return logical 
#' @export
#' @examples
is_tip <- function(node, tree){
  if (node==length(tree$tip.label)+1){ # check if node is root 
    bool <- FALSE
  }
  else{
    df <- tree$orig
    bool <- !df[df["idx2"]==node,]$split # has the node children?
  }
  return(bool)
}


#' Get the two children of a node
#'
#' The left children is always the node the further from the root
#' The right children is always the node the closest to the root
#' If the node is a tip returns NA values 
#'
#' @param node integer, index of the node
#' @param tree phylo tree
#'
#' @return list, $left = left child, $right = right child
#' @export
#' @examples
get_child <- function(node, tree){
  child <- list("left" = NA, "right" = NA)
  dist.all <- castor::get_all_distances_to_root(tree)
  
  # First check that the node is not a tip (otherwise returns NA)
  if (!is_tip(node, tree)){ 
    df  <- tree$orig
    idx <- df[df["parent2"]==node,]$idx2
    
    # If both child are tips (same dist. to root), random attribution (l/r)
    if (is_tip(idx[1], tree) & is_tip(idx[2], tree)){
      child$left  <- idx[1]
      child$right <- idx[2]
    }
    
    # Else (>= 1 child is internal), left child is the further from root
    else {
      dist.child <- c(dist.all[idx[1]], dist.all[idx[2]])
      left  <- which(dist.child == max(dist.child))
      right <- which(dist.child == min(dist.child))
      child$left  <- idx[left]
      child$right <- idx[right]
    }
  }
  return(child) 
}


#' Traverse a tree inorder 
#' 
#'
#' @param tree phylo tree
#'
#' @return vector, ordered sequence of node indexes of the inorder traversal^
#' @export
#' @examples
traverse_inorder <- function(tree){
  node  <- length(tree$tip.label) + 1 # root index 
  stack <- c()
  inorder <- c()
  while(length(stack) != 0 | !is.na(node)){
    if(!is.na(node)){
      stack <- c(node, stack)
      node <- get_child(node, tree)$left
    }
    else{
      node <- stack[1]
      inorder <- c(inorder, node)
      stack <- stack[stack!=node]
      node <- get_child(node, tree)$right
    }
  }
  return(inorder)
}


#' Compute all the distances of tips to their most recent ancestor 
#'
#' For each node, 
#' 1. find its parent (most recent ancestor)
#' 2. compute the distance between the node and its parent 
#' This distance is only computed for tips, as we don't need this distance 
#' for the internal nodes in the encoding. 
#' Thus internal node indexes are filled with NA 
#' See Voznica 2021 - bioRxiv - DOI:10.1101/2021.03.11.435006 
#'
#' @param tree phylo tree
#'
#' @return vector of distance 
#' @export
#' @examples
get_all_distances_to_ancestor <- function(tree){
  n <- length(tree$orig$idx) + 1 
  n_node <- tree$Nnode
  n_tip <- n - n_node
  dist <- c()
  
  # Fill tips value w/ their distance to their parent node
  for (i in 1:n_tip){
    df <- tree$orig
    df.i <- df[df["idx2"]==i,]
    parent <- df.i$parent2
    bool <- tree$edge[,1]==parent & tree$edge[,2]==i
    edge.idx <- which(bool==TRUE)
    edge.length <- tree$edge.length[edge.idx]
    dist <- c(dist, edge.length)
  }
  
  # Fill nodes values w/ NA
  for (i in 1:n_node){ 
    dist <- c(dist, NA)
  }
  return(dist)
}

#' Encode a phylogenetic tree into a vector 
#'
#' Compute the compact encoding of a phylogenetic tree. This encoding is
#' bijective. The encoding methods is named: "Compact Bijective Ladderized 
#' Vector" (CBLV).
#' See Voznica 2021 - bioRxiv - DOI:10.1101/2021.03.11.435006
#'
#' @param tree phylo tree to encode 
#'
#' @return vector containing the encoding
#' @export
#' @examples
encode_phylo <- function(tree){
  # 
  # A list for nodes containing their distances to the root 
  # A list for tips containing their distances to the most recent ancestor 
  inorder <- traverse_inorder(tree)
  tips  <- c()
  nodes <- c()
  dist_to_root <- castor::get_all_distances_to_root(tree)
  dist_to_ancestor <- get_all_distances_to_ancestor(tree)
  for (node in inorder){
    if (is_tip(node, tree)){
      tips <- c(tips, dist_to_ancestor[node])
    }
    else{
      nodes <- c(nodes, dist_to_root[node])
    }
  }
  encoding <- list("nodes" = nodes, "tips" = tips)
  return(encoding)
}


#### end ####

if (FALSE){
  tree <- trees(c(.1, 0), "bd", max.taxa=10)[[1]]
  plot(tree, edge.width = 2, label.offset = 0.1, type = "cladogram")
  nodelabels()
  tiplabels()
  
  traverse_inorder(tree)
  encode_phylo(tree)
}


