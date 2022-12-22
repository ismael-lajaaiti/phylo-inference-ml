# Libraries and Sources
source("R/libraries.R")
source("R/phylogeny-properties.R")
source("R/phylogeny-graph-data.R")
source("infer-general-functions.R")


n_core <- 10
n_rep <- 50
prefix_raw_phylo <- "bisse-1e6-phylogenies/raw/bisse-n20000-trees-and-params"
df_list <- mclapply(1:n_rep, function(i) {
    out <- readRDS(paste(prefix_raw_phylo,
        str_pad(i, 2, pad = "0"), ".rds",
        sep = ""
    ))
    tree_list <- out$trees
    n_tree <- 10 # length(tree_list)
    list_tree_df <- list()
    for (k in 1:n_tree) {
        tree <- tree_list[[k]]
        edge <- edge_df(tree)
        node <- node_df(tree)
        df <- list("node" = node, "edge" = edge)
        list_tree_df[[k]] <- df
    }
}, mc.cores = n_core)

# Saving
dir <- "trees-dataset"
fname <- get_backbone_save_name(n_tree, n_taxa, param.range) # save name for file
fname <- paste(fname, "sscheck", ss_check, "df.rds", sep = "-")
fname <- paste(dir, fname, sep = "/")
saveRDS(list.df.tree, fname) # saving file
cat(paste(fname, " saved.\n"))
