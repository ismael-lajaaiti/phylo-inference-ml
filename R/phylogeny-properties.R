#' Is node the root of the tree?
is_root <- function(tree, node_idx) {
    node_idx == ape::Ntip(tree) + 1
}

#' Is node a tip of the tree?
is_tip <- function(tree, node_idx) {
    node_idx <= ape::Ntip(tree) & node_idx >= 1
}
