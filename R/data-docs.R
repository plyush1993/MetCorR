#' Example intensity matrix
#'
#' A small LC–MS intensity table for examples and tests.
#'
#' @format A data frame with N rows (samples) and P numeric feature columns.
#'   Row names are sample IDs.
#' @examples
#' data(example_intensity, package = "MetCorR")
#' example_intensity[c(1:5), c(1:2)]
"example_intensity"

#' Example sample metadata
#'
#' @format A metadata data frame associated with example_intensity with columns:
#' \describe{
#'   \item{class}{character}
#'   \item{order}{integer}
#'   \item{batch}{integer}
#' }
#' @examples
#' data(example_meta, package = "MetCorR")
#' head(example_meta)
"example_meta"
