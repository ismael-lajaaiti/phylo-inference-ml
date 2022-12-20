# Convert phylogenies from the tree (ape) format to LTT and save.

source("R/libraries.R")
source("infer-general-functions.R")
source("R/neural-network-functions.R")

n_core <- 10
n_rep <- 50
range_tree_size <- c(100, 1000)
prefix_raw_phylo <- "bisse-1e6-phylogenies/raw/bisse-n20000-trees-and-params"
df_list <- mclapply(1:n_rep, function(i) {
    out <- readRDS(paste(prefix_raw_phylo,
        str_pad(i, 2, pad = "0"), ".rds",
        sep = ""
    ))
    tree_list <- out$trees
    generate_ltt_dataframe(
        tree_list,
        range_tree_size,
        verbose = FALSE
    )
}, mc.cores = n_core)
df_ltt <- do.call(cbind, df_list)

saveRDS(df_ltt, "bisse-1e6-phylogenies/formatted/ltt/bisse-ltt-all.rds")
