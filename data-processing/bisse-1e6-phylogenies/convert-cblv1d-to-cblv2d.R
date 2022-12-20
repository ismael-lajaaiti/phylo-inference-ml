# Convert and save the phylogenies from 1D-CBLV encoding to the 2D-CBLV.
# ----------------------------------------------------------------------

# Load dependencies.
source("infer-general-functions.R")

# Load the CBLV 1d matrix.
fname_cblv_1d <- "bisse-1e6-phylogenies/formatted/cblv/bisse-cblv-all.rds"
cblv_1d <- readRDS(fname_cblv_1d)

# From 1d to 2d.
flat_vec <- c(cblv_1d) # (3_000,1_000_000) -> (3_000_000_000,)
cblv_2d <- matrix(flat_vec, 3 * ncol(enc_1d), 1000, byrow = TRUE)
fname_cblv_2d <- "bisse-1e6-phylogenies/formatted/cblv/bisse-cblv2d-all.rds"

saveRDS(cblv_2d, fname_cblv_2d)
