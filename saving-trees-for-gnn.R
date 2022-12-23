source("R/libraries.R")
source("R/phylogeny-properties.R")
source("R/phylogeny-graph-data.R")
source("infer-general-functions.R")

n_core <- 10
n_rep <- 50
prefix_raw_phylo <- "bisse-1e6-phylogenies/raw/bisse-n20000-trees-and-params"
prefix_gnn <- "bisse-1e6-phylogenies/formatted/graph/bisse-n20000-node-edge"
mclapply(1:n_rep, function(i) {
    out <- readRDS(
        stringr::str_c(
            prefix_raw_phylo,
            stringr::str_pad(i, 2, pad = "0"),
            ".rds"
        )
    )
    tree_list <- out$trees
    n_tree <- length(tree_list)
    df_list <- list()
    for (k in 1:n_tree) {
        tree <- tree_list[[k]]
        edge <- edge_df(tree)
        node <- node_df(tree)
        df <- list("node" = node, "edge" = edge)
        df_list[[k]] <- df
    }
    fname <- stringr::str_c(
        prefix_gnn,
        stringr::str_pad(i, 2, pad = "0"),
        ".rds"
    )
    saveRDS(df_list, fname)
    print(str_c(fname, " saved."))
}, mc.cores = n_core)
