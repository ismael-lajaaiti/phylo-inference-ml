# Remove tips from the 1D-CBLV encoding and save.
# -----------------------------------------------

source("infer-general-functions.R")

"Remove tips in the 1d CBLV encoding."
remove_tip_state <- function(cblv_matrix, verbose = FALSE) {
    cblv_matrix[1:2000, ]
}

cblv_directory <- "bisse-1e6-phylogenies/formatted/cblv/"
cblv_matrix <- readRDS(str_c(cblv_directory, "bisse-cblv-all.rds"))
cblv_matrix_notipstate <- remove_tip_state(cblv_matrix)

saveRDS(
    cblv_matrix_notipstate,
    str_c(cblv_directory, "bisse-cblv-notipstate-all.rds")
)
