# This file contains the functions to compute the summary statistics 
# of a phylogenetic tree
# See Saulnier Et. Al. 2017 - PLOS Computational Biology
# DOI:10.1371/journal.pcbi.1005416


#### Libraries & Sources####


library(ape)
library(diversitree)
library(castor)
library(phangorn)
source("encode-cblv.R")


#### end ####


#### Nodes data.frame ####


create_nodes_data_frame <- function(tree){
  # Initialize the data.frame containing nodes variables 
  n <- length(tree$orig$idx) + 1 
  n_tips <- length(tree$tip.label)
  df.nodes <- data.frame(index=rep(NA,n),   # index of the node
                         is.tip=rep(NA,n),  # is the node a tip?
                         dist=rep(NA,n),    # distance to the root
                         part=rep(NA,n),    # (3*dist) %/% h + 1 in {1;2;3}
                         depth=rep(NA,n),   # distance the root (in edges)
                         colless=rep(NA,n), # colless score
                         stair=rep(NA,n))   # staircaseness 
  return(df.nodes)
}


fill_nodes_index <- function(df.nodes){
  # Fill the "index" column of the nodes data.frame
  n <- length(df.nodes[,1])
  df.nodes["index"]  <- 1:n
  return(df.nodes)
}


fill_nodes_is.tip <- function(df.nodes, tree){
  # Fill the "is.tip" column of the nodes data.frame
  n_tips <- length(tree$tip.label)
  df.nodes["is.tip"] <- as.logical(df.nodes["index"] <= n_tips)
  return(df.nodes)
}


fill_nodes_dist <- function(df.nodes){
  # Fill the "dist" column of the nodes data.frame
  df.nodes["dist"]   <- castor::get_all_distances_to_root(tree)
  return(df.nodes)
}


fill_nodes_part <- function(df.nodes){
  # Fill the "part" column of the nodes data.frame
  dist <- df.nodes$dist 
  height <- get_height(df.nodes)
  greater1 <- as.numeric(dist > height/3)
  greater2 <- as.numeric(dist > as.numeric(2*height/3))
  df.nodes["part"] <- 1 + greater1 + greater2 
  return(df.nodes)
}


fill_nodes_colless <- function(df.nodes, tree){
  # Fill the "colless" column of the nodes data.frame
  n <- length(df.nodes[,1])
  for (i in 1:n){
    if (df.nodes[i, "is.tip"] == FALSE){
      df.nodes[i, "colless"] <- get_nodes_colless(df.nodes, tree, i)
    }
  }
  return(df.nodes)
}


fill_nodes_stair <- function(df.nodes, tree){
  # Fill the "stair" column of the nodes data.frame
  n <- length(df.nodes[,1])
  for (i in 1:n){
    if (df.nodes[i, "is.tip"] == FALSE){
      df.nodes[i, "stair"] <- get_nodes_staircaseness2(df.nodes, tree, i)
    }
  }
  return(df.nodes)
}


fill_nodes_depth <- function(df.nodes, tree){
  # Fill the "depth" column of the nodes data.frame
  dist_edge <- castor::get_all_distances_to_root(tree, as_edge_count=TRUE)
  df.nodes["depth"] <- dist_edge
  return(df.nodes)
}


fill_nodes_all <- function(df.nodes, tree){
  # Fill all the columns of the nodes data.frame 
  df.nodes <- fill_nodes_index(df.nodes)
  df.nodes <- fill_nodes_is.tip(df.nodes, tree)
  df.nodes <- fill_nodes_dist(df.nodes)
  df.nodes <- fill_nodes_part(df.nodes)
  df.nodes <- fill_nodes_depth(df.nodes, tree)
  df.nodes <- fill_nodes_colless(df.nodes, tree)
  df.nodes <- fill_nodes_stair(df.nodes, tree)
  return(df.nodes)
}


get_nodes_colless <- function(df.nodes, tree, node){
  # Get the colless score of a given node of a tree
  
  n_child_list <- list(0, 0) # store the number of left & right children
  children <- phangorn::Children(tree, node)
  
  for (i in 1:2){
    desc_child <- phangorn::Descendants(tree, children[i])[[1]]
    n_child <- length(desc_child) # number of children on one side (L or R)
    n_child_list[[i]] <- n_child  # save the value 
  }
  
  # colless = |number of left children - number of right children|
  colless <- abs(n_child_list[[1]] - n_child_list[[2]])
  
  return(colless)
}


get_nodes_staircaseness2 <- function(df.nodes, tree, node){
  # Get the stairecaseness of a given of a tree
  
  n_child_list <- c(0, 0) # store the number of left & right children
  children <- phangorn::Children(tree, node)
  
  for (i in 1:2){
    desc_child <- phangorn::Descendants(tree, children[i])[[1]]
    n_child <- length(desc_child) # number of children on one side (L or R)
    n_child_list[[i]] <- n_child  # save the value 
  }
  
  # stair = (min(left or right children)/max(left or right children))
  stair <- min(n_child_list)/max(n_child_list)
  
  return(stair)
}

#### end ####


#### Edges data.frame ####


create_edges_data_frame <- function(tree){
  # Initialize the data.frame containing nodes variables 
  n <- length(tree$edge.length)
  df.edges <- data.frame(index=rep(NA,n),   # index of the edge  
                         is.ext=rep(NA,n),  # is the branch external? 
                         length=rep(NA,n),  # length of the branch
                         part=rep(NA,n),    # part of the tree (n/3)
                         node1=rep(NA,n),   # first node (older)
                         node2=rep(NA,n))   # second node (younger)
  return(df.edges)
}


fill_edges_index <- function(df.edges){
  # Fill the "index" column of the edges data.frame
  n <- length(df.edges[,1])
  df.edges["index"]  <- 1:n
  return(df.edges)
}


fill_edges_is.ext <- function(df.edges, df.nodes){
  # Fill the "is.ext" column of the edges data.frame
  df.edges["is.ext"] <- df.nodes[df.edges$node2, "is.tip"]
  return(df.edges)
}


fill_edges_length <- function(df.edges){
  # Fill the "length" column of the edges data.frame
  df.edges["length"] <- tree$edge.length
  return(df.edges)
}


fill_edges_part <- function(df.edges, df.nodes){
  # Fill the "part" column of the edges data.frame
  df.edges["part"] <- df.nodes[df.edges$node2, "part"]
  return(df.edges)
}


fill_edges_node12 <- function(df.edges){
  # Fill the "node1" & "node2" column of the edges data.frame
  df.edges[c("node1","node2")] <- tree$edge
  return(df.edges)
}


fill_edges_all <- function(df.edges, df.nodes){
  # Fill all the columns of the edges data.frame
  df.edges <- fill_edges_index(df.edges)
  df.edges <- fill_edges_length(df.edges)
  df.edges <- fill_edges_node12(df.edges)
  df.edges <- fill_edges_part(df.edges, df.nodes)
  df.edges <- fill_edges_is.ext(df.edges, df.nodes)
  return(df.edges)
}


#### end ####


#### Summary Statistics  ####


# bl := branch length 
# a  := all branches 
# e  := external branches (i.e. branches connected to a tip)
# i  := piecewise branches i.e. branches in the part {1/3;2/3;3/3}
# ie := ratio of piecewise internal branches over external branches
# stats stands for mean, median and variance 


get_bl_stats <- function(bl){
  # Returns the mean, median and variance of a list of branch lengths (bl)
  stats <- list("mean"=-999, "med"=-999, "var"=-999)
  stats$mean <- mean(bl)
  stats$med  <- median(bl)
  stats$var  <- var(bl)
  return(stats)
}


get_bl_ss <- function(df.edges){
  # Returns the list of all branch lengths summary statistics (24)
  
  ss.bl   <-   list("a.mean"   =-999, "a.med"   =-999, "a.var"   =-999,
                    "e.mean"   =-999, "e.med"   =-999, "e.var"   =-999,
                    "i.1.mean" =-999, "i.1.med" =-999, "i.1.var" =-999,
                    "i.2.mean" =-999, "i.2.med" =-999, "i.2.var" =-999,
                    "i.3.mean" =-999, "i.3.med" =-999, "i.3.var" =-999,
                    "ie.1.mean"=-999, "ie.1.med"=-999, "ie.1.var"=-999,
                    "ie.2.mean"=-999, "ie.2.med"=-999, "ie.2.var"=-999,
                    "ie.3.mean"=-999, "ie.3.med"=-999, "ie.3.var"=-999)
  
  df.edges.ext <- df.edges[df.edges["is.ext"]==TRUE,]
  df.edges.int <- df.edges[df.edges["is.ext"]==FALSE,]
  
  bl.a <- df.edges$length
  bl.e <- df.edges.ext$length
  bl.i.1 <- df.edges.int[df.edges.int["part"]==1,]$length
  bl.i.2 <- df.edges.int[df.edges.int["part"]==2,]$length
  bl.i.3 <- df.edges.int[df.edges.int["part"]==3,]$length
  
  bl.list <- list(bl.a, bl.e, bl.i.1, bl.i.2, bl.i.3)
  n  <- length(bl.list)
  
  for (i in 1:n){
    bl <- bl.list[[i]]
    stats <- get_bl_stats(bl)
    ss.bl[[3*(i-1)+1]] <- stats$mean
    ss.bl[[3*(i-1)+2]] <- stats$med
    ss.bl[[3*(i-1)+3]] <- stats$var
  }
  
  for (i in 1:3){
    for (stat in c("mean", "med", "var")){
      name.i  <- paste("i" , i, stat, sep=".")
      name.ie <- paste("ie", i, stat, sep=".")
      name.e <- paste("e", stat, sep=".")
      stat.e <- ss.bl[[name.e]]
      stat.i <- ss.bl[[name.i]]
      ratio <- stat.i / stat.e
      ss.bl[[name.ie]] <- ratio
    }
  }
  return(ss.bl)
}


is_in_ladder <- function(node, tree){
  sum <- 0
  if (!is_tip(node, tree)){
    children <- phangorn::Children(tree, node)
    sum <- as.numeric(is_tip(children[1], tree) + is_tip(children[2], tree))  
  }
  return(sum==1)
}


traverse_ladder <- function(tree, node, preorder){
  n_tips <- length(tree$tip.label)
  root <- n_tips + 1 
  if (is_in_ladder(node, tree)){
    preorder <- c(preorder, node)
  }
  else{
    preorder <- c(preorder, -1)
  }
  children <- get_child(node, tree)
  child.left <- children$right
  child.right <- children$left
  if (!is.na(child.left)){
    preorder <- traverse_ladder(tree, child.left, preorder)
  }
  if (!is.na(child.right)){
    preorder <- traverse_ladder(tree, child.right, preorder)
  }
  return(preorder)
}


extract_ladders <- function(traverse_ladder){
  n <- length(traverse_ladder)
  ladder.list <- list()
  ladder  <- c()
  for (i in 1:n){
    x <- traverse_ladder[i]
    if (x != -1){ladder <- c(ladder, x)}
    else if (x==-1 & length(ladder)==1){ladder <- c()}
    else if (x==-1 & length(ladder)>=2){
      ladder.list <- append(ladder.list, list(ladder))
      ladder <- c()
    }
  }
  return(ladder.list)
}


get_ladders <- function(tree){
  root <- length(tree$tip.label) + 1 
  traverse_ladder <- traverse_ladder(tree, root, c())
  ladders <- extract_ladders(traverse_ladder)
  return(ladders)
}


get_max_ladder <- function(tree){
  max_l <- 0
  n_tips <- length(tree$tip.label)
  ladders <- get_ladders(tree)
  for (ladder in ladders){
    l <- length(ladder)
    if (l > max_l){max_l <- l}
  }
  return(max_l / n_tips)
}


get_il_nodes <- function(tree){
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips - 1 
  ladders <- get_ladders(tree)
  il_nodes <- 0
  for (ladder in ladders){
    l <- length(ladder)
    il_nodes <- il_nodes + l 
  }
  return(il_nodes / n_nodes)
}


get_staircaseness1 <- function(df.nodes){
  # Returns the proportion of imbalanced internal nodes 
  # that have different numbers of leaves between the left and the right side
  n <- length(tree$orig$idx) + 1 
  n_tips <- length(tree$tip.label)
  n_nodes <- n - n_tips
  stair1 <- table(df.nodes[df.nodes["stair"]==1,]$stair)[[1]] / n_nodes 
  return(1 - stair1)
}


get_topo_ss <- function(df.nodes){
  # Returns the list of all topological summary statistics excluding bl (6)
  
  ss.topo <- list("height"  = -999, "colless" = -999, "sackin"     = -999, 
                  "wdratio" = -999, "deltaw"  = -999, "max_ladder" = -999,
                  "il_nodes"= -999, "stair1"  = -999, "stair2"     = -999)
  
  ss.topo$height     <- get_height(df.nodes)
  ss.topo$colless    <- df.nodes[df.nodes["is.tip"]==FALSE,]$colless
  ss.topo$sackin     <- df.nodes[df.nodes["is.tip"]==TRUE,]$depth
  ss.topo$wdratio    <- get_wdratio(df.nodes) 
  ss.topo$deltaw     <- get_deltaw(df.nodes)
  ss.topo$max_ladder <- get_max_ladder(tree)
  ss.topo$il_nodes   <- get_il_nodes(tree)
  ss.topo$stair1     <- get_staircaseness1(df.nodes) 
  ss.topo$stair2     <- df.nodes[df.nodes["is.tip"]==FALSE,]$stair
  return(ss.topo)
}


get_height <- function(df.nodes){
  # Returns the height of the tree (i.e. max distance to root)
  height <- max(df.nodes$dist)
  return(height)
}


get_wdratio <- function(df.nodes){
  # Returns the ratio of the maximal width over the maximal depth 
  width <- max(table(df.nodes$depth))
  depth <- max(df.nodes$depth)
  ratio <- width/depth
  return(ratio)
}


get_deltaw <- function(df.nodes){
  # Returns the maximal difference in width between two consecutive depths
  table.width <- table(df.nodes$depth)
  n <- length(table.width)-1
  diff <- 0
  for (i in 1:n){
    diff <- max(diff, abs(table.width[[i]] - table.width[[i+1]]))
  }
  return(diff)
}


get_ltt_coords <- function(tree){
  # Returns 20 coordinates (x=time, y=number of lineages) of the LTT 
  ltt.coord <- ape::ltt.plot.coords(tree) 
  n_events <- length(ltt.coord[,1])
  step <- n_events %/% 20
  bins <- c()
  for (i in 0:19){
    bins <- c(bins, 1 + i*step)
  }
  out.ltt.coord <- ltt.coord[bins,]
  return(out.ltt.coord)
}


get_ltt_slopes <- function(tree){
  
  slopes <- list("p.1" = -999, "p.2" = -999, "p.3" = -999, 
                 "ratio.12" = -999, "ratio.23" = -999)
  
  ltt.coord <- ape::ltt.plot.coords(tree) 
  ltt.df <- as.data.frame(ltt.coord)
  ltt.df["N"] <- log(ltt.df["N"])
  t <- ltt.df$time[1] # root time 
  ltt.df.1 <- ltt.df[ltt.df["time"] < 2*t/3,]
  ltt.df.2 <- ltt.df[ltt.df["time"] <= t/3 & ltt.df["time"] >= 2*t/3,]
  ltt.df.3 <- ltt.df[ltt.df["time"] > t/3,]
  slopes$p.1  <- lm(N ~ time, ltt.df.1)$coef[[2]]
  slopes$p.2  <- lm(N ~ time, ltt.df.2)$coef[[2]]
  slopes$p.3  <- lm(N ~ time, ltt.df.3)$coef[[2]]
  slopes$ratio.12 <- slopes$p.2 / slopes$p.1
  slopes$ratio.23 <- slopes$p.3 / slopes$p.2
  return(slopes)
}


get_all_ss <- function(df.nodes, df.edges, tree){
  # Returns a list with all the summary statistics (ss) of the tree
  
  ss.bl   <- get_bl_ss(df.edges)                 # branch length ss 
  ss.topo <- get_topo_ss(df.nodes)               # topological ss
  ss.ltt <- list("coord" = -999, slopes = -999)  # ltt ss
  ss.ltt$coord <- get_ltt_coords(tree)
  ss.ltt$slopes <- get_ltt_slopes(tree)
  
  ss.all <- c(ss.bl, ss.topo, ss.ltt)    # merge all ss in a single list 
  return(ss.all)
}


get_ss <- function(tree){
  # Returns a list of all the summary statistics, given a tree directly
  
  # Create nodes data.frame
  df.nodes <- create_nodes_data_frame(tree)
  df.nodes <- fill_nodes_all(df.nodes, tree)
  
  # Create edges data.frame
  df.edges <- create_edges_data_frame(tree)
  df.edges <- fill_edges_all(df.edges, df.nodes)
  
  # Compute all summary statistics 
  ss <- get_all_ss(df.nodes, df.edges, tree)
  return(ss)
}


#### end ####


ggplot(ltt.df.3, aes(x=time, y=N))+
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  stat_regline_equation()


tree <- trees(c(.1, 0), "bd", max.taxa=30)[[1]]
plot(tree)
nodelabels()
tiplabels()

root <- length(tree$tip.label) + 1 
traverse_ladder <- traverse_ladder(tree, root, c())
ladders <- get_ladders(traverse_ladder)
for (l in ladders){
  print(length(l))
}



df.nodes

get_ss(tree)

df.nodes <- create_nodes_data_frame(tree)
df.nodes <- fill_nodes_all(df.nodes, tree)
df.edges <- create_edges_data_frame(tree)
df.edges <- fill_edges_all(df.edges, df.nodes)
