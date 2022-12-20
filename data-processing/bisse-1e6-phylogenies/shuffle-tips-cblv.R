# Shuffle tips from the 1D-CBLV encoding and save.
# ------------------------------------------------

source("infer-general-functions.R")

"Shuffle tip states in the 1d CBLV encoding."
shuffle_tip_state <- function(cblv_matrix, ntips_trees, verbose = FALSE) {
    n_trees <- length(ntips_trees)
    for (i in 1:n_trees) {
        n_tips <- ntips_trees[i]
        end <- 2000 + n_tips
        cblv_matrix[2000:end, i] <- sample(cblv_matrix[2000:end, i])
        if (verbose) { # print progression bar
            svMisc::progress(i,
                max.value = n_trees,
                progress.bar = TRUE, init = (i == 1)
            )
        }
    }
    cblv_matrix
}

cblv_directory <- "bisse-1e6-phylogenies/formatted/cblv/"
cblv_matrix <- readRDS(str_c(cblv_directory, "bisse-cblv-all.rds"))
ntips_trees <- readRDS("bisse-1e6-phylogenies/raw/ntips-all.rds")
cblv_matrix_shuffled <- shuffle_tip_state(cblv_matrix, ntips_trees)

saveRDS(
    cblv_matrix_shuffled,
    str_c(cblv_directory, "bisse-cblv-shuffled-all.rds")
)
