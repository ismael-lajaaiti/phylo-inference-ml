#' Extract given indices in a list of vectors.
extract_elements <- function(list_of_vecs, indices) {
    as.list(do.call(cbind, list_of_vecs)[indices, ] |> as.data.frame())
}

#' Add two columns to the summary statistic dataframe corresponding to
#' the diversification parameteres.
add_params_to_sumstat <- function(df, params) {
    for (name in names(params)) {
        df[, name] <- params[[name]]
    }
    df
}
