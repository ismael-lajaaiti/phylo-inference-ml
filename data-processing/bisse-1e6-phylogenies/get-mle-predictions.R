# Extract MLE predictions on the phylogenies of the test set.

source("R/libraries.R")
source("R/utils.R")

dir <- "bisse-1e6-phylogenies/"
mle_predictions <- list(lambda = c(), q = c())
for (i in 1:50) {
    temp <- readRDS(
        str_c(
            dir,
            "mle/bisse-n20000-predmle",
            str_pad(i, 2, pad = "0"),
            ".rds"
        )
    )
    mle_predictions$lambda <- c(mle_predictions$lambda, temp$lambda0)
    mle_predictions$q <- c(mle_predictions$q, temp$q01)
}

test_indices <- readRDS(str_c(dir, "indices-1e6-phylogenies.rds"))$test
mle_predictions_test <- extract_elements(mle_predictions, test_indices)
saveRDS(mle_predictions, str_c(dir, "mle/bisse-predmle-all.rds"))
saveRDS(mle_predictions_test, str_c(dir, "predictions/mle-predictions.rds"))
