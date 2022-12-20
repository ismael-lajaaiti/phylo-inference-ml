# Compute the number of tips of each phylogeny and save.
# ------------------------------------------------------

ntips_trees <- c() # initilization

number_of_tips <- function(tree) tree$Nnode + 1
prefix_fname <- "bisse-1e6-phylogenies/raw/bisse-n20000-trees-and-params"

for (i in 1:50) {
    print(i) # print progress
    out <- readRDS(str_c(prefix_fname, str_pad(i, 2, pad = "0"), ".rds"))
    trees <- out$trees
    ntips_trees_temp <- sapply(trees, number_of_tips)
    ntips_trees <- c(ntips_trees, ntips_trees_temp)
}

saveRDS(ntips_trees, "bisse-1e6-phylogenies/raw/ntips-all.rds")
