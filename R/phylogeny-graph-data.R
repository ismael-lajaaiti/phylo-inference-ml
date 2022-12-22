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
edge_df <- function(tree) {
    df <- as.data.frame(tree$edge)
    colnames(df) <- c("node1", "node2")
    df
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
node_df <- function(tree) {
    n_nodes <- 2 * tree$Nnode + 1
    attr_name <- c(
        "ancestor",
        "clade_asym",
        "descendant",
        "dist",
        "mean_edge",
        "state",
        "time_asym"
    )
    n_attr <- length(attr_name)
    attr_nodes <- vector(mode = "list", length = n_attr)
    names(attr_nodes) <- attr_name
    attr_nodes$dist <- castor::get_all_distances_to_root(tree,
        as_edge_count = FALSE
    )
    attr_nodes$ancestor <- castor::get_all_distances_to_root(tree,
        as_edge_count = TRUE
    )
    for (i in 1:n_nodes) {
        descendant <- length(phangorn::Descendants(tree, i, "all")) - 1
        mean_edge <- mean_edge_length(tree, i)
        time_asym <- time_asym(tree, i)
        clade_asym <- clade_asym(tree, i)
        state <- node_state(tree, i)
        attr_nodes$descendant <- c(attr_nodes$descendant, descendant)
        attr_nodes$clade_asym <- c(attr_nodes$clade_asym, clade_asym)
        attr_nodes$time_asym <- c(attr_nodes$time_asym, time_asym)
        attr_nodes$mean_edge <- c(attr_nodes$mean_edge, mean_edge)
        attr_nodes$state <- c(attr_nodes$state, state)
    }
    as.data.frame(attr_nodes)
}

in_edge_length <- function(tree, node) {
    if (is_root(tree, node)) {
        return(0.0)
    }
    parent <- phangorn::Ancestors(tree, node, "parent") # get node's parent
    bool <- tree$edge[, 1] == parent & tree$edge[, 2] == node
    index <- which(bool == TRUE) # edge index
    tree$edge.length[index] # edge length
}

mean_edge_length <- function(tree, node) {
    # Node is root - mean of two outcoming edges
    if (is_root(tree, node)) {
        children <- phangorn::Children(tree, node) # get children
        length1 <- in_edge_length(tree, children[1]) # length: root --> child1
        length2 <- in_edge_length(tree, children[2]) # length: root --> child2
        mean <- mean(c(length1, length2)) # compute the mean
    } else if (is_tip(tree, node)) {
        mean <- in_edge_length(tree, node)
    } else {
        children <- phangorn::Children(tree, node) # get children
        length1 <- in_edge_length(tree, children[1]) # length: node --> child1
        length2 <- in_edge_length(tree, children[2]) # length: node --> child2
        length3 <- in_edge_length(tree, node) # length: parent --> node
        mean <- mean(c(length1, length2, length3)) # compute the mean
    }
    mean
}

clade_asym <- function(tree, node) {
    if (is_tip(tree, node)) {
        return(0.0)
    }
    children <- phangorn::Children(tree, node) # get children
    length1 <- in_edge_length(tree, children[1]) # length: node --> child1
    length2 <- in_edge_length(tree, children[2]) # length: node --> child2
    abs(length1 - length2) / mean(c(length1, length2))
}

time_asym <- function(tree, node) {
    if (is_tip(tree, node) || is_root(tree, node)) {
        return(0.0)
    }
    l1 <- in_edge_length(tree, node) # length: parent --> node
    children <- phangorn::Children(tree, node) # get children
    l2 <- in_edge_length(tree, children[1]) # length: node --> child1
    l3 <- in_edge_length(tree, children[2]) # length: node --> child2
    abs(l1 - mean(c(l2, l3))) / mean(c(l1, l2, l3))
}

node_state <- function(tree, node) {
    if (is_tip(tree, node)) { # node is tip, save state
        state <- tree$tip.state[[node]]
    } else { # node isn't tip, state is unknown (i.e. -1)
        state <- -1
    }
    state
}
