# Encode phylogeny with the CBLV-encoding and save.
# -------------------------------------------------

# Load mandatory dependencies
source("infer-general-functions.R")
source("neural-network-functions.R")

n_taxa <- c(100, 1000) # phylogeny range size

# Define utility functions.
"Replace `old_val` by `new_val` in a given vector `vec`."
replace_values <- function(vec, old_val, new_val) {
    n <- length(vec)
    for (i in 1:n) {
        if (vec[i] == old_val) {
            vec[i] <- new_val
        }
    }
    return(vec)
}
"Convert tip states as follow: 0 -> -1, 1 -> 1."
relabel_tips <- function(mat, trees) {
    for (i in 1:length(trees)) {
        n_tips <- trees[[i]]$Nnode + 1
        mat[2000:(2000 + n_tips), i] <- replace_values(
            mat[2000:(2000 + n_tips), i],
            0, -1
        )
    }
    return(mat)
}

# Read raw phylogenies.
prefix_raw_phylo <- "bisse-1e6-phylogenies/raw/bisse-n20000-trees-and-params"
prefix_cblv <- "bisse-1e6-phylogenies/formatted/cblv/bisse-n20000-cblv"
out <- readRDS(paste(prefix_raw_phylo,
    str_pad(1, 2, pad = "0"), ".rds",
    sep = ""
))
trees <- out$trees # extract phylogenies

# Read and process CBLV encoding.
cblv_matrix <- readRDS(paste(prefix_cblv,
    str_pad(1, 2, pad = "0"), ".rds",
    sep = ""
))
cblv_matrix <- relabel_tips(cblv_matrix, trees)

# Repeat on data files.
for (i in 2:50) {
    print(i)
    out <- readRDS(paste(prefix_raw_phylo,
        str_pad(i, 2, pad = "0"), ".rds",
        sep = ""
    ))
    trees <- out$trees
    cblv_matrix_temp <- readRDS(paste(prefix_cblv,
        str_pad(i, 2, pad = "0"), ".rds",
        sep = ""
    ))
    cblv_matrix_temp <- relabel_tips(cblv_matrix_temp, trees)
    cblv_matrix <- cbind(cblv_matrix, cblv_matrix_temp)
}

saveRDS(cblv_matrix, "bisse-1e6-phylogenies/formatted/cblv/bisse-cblv-all.rds")
