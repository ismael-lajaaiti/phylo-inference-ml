#' Extract given indices in a list of vectors.
extract_elements <- function(list_of_vecs, indices) {
    as.list(do.call(cbind, list_of_vecs)[indices, ] |> as.data.frame())
}
