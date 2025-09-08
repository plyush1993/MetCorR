#' MetCorR: GAM-based batch/order correction for metabolomics peak table
#'
#' @description
#' Corrects intensity matrices using GAM fits on QC samples.
#' Two modes:
#' 1) y ~ s(order)
#' 2) y ~ s(order, batch)
#'
#' @param method Integer; 1 or 2. Default 2. 1 - one feature mode, only run order; 2 - two features mode, run order & batch index.
#' @param int_data A data.frame of intensities (rows = samples, cols = features).
#' @param order Numeric or coercible; injection order per sample (length = nrow(int_data)).
#' @param class  Character or factor; sample class/labels per sample (length = nrow(int_data)).
#' @param batch  Numeric or coercible; batch index per sample (length = nrow(int_data)). Required for method 2.
#' @param qc_label Pattern used to identify QC rows in class via \code{grep(qc_label, class)}.
#'
#' @return A numeric data.frame of corrected intensities (same dimnames as input).
#' @examples
#' \dontrun{
#'  data(example_intensity, package = "MetCorR")
#'  data(example_meta, package = "MetCorR")
#'  out <- MetCorR(
#'    method   = 2,
#'    int_data = example_intensity,
#'    order    = example_meta$order,
#'    class    = example_meta$class,
#'    batch    = example_meta$batch,
#'    qc_label = "QC"
#'  )
#' }
#'
#' @importFrom mgcv gam s
#' @importFrom pbapply pblapply pboptions
#' @importFrom crayon blue red cyan green %+%
#' @export
MetCorR <- function(method = 2, int_data, order, class, batch, qc_label) {

  # Banner
  cat("\n")
  cat(crayon::red(
    c("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
  ),
  sep = "\n")

  cat(crayon::blue(
    c(
      "  .--.   ,--.         ,--.   ,-----.              ,------.    ",
      "  |   '.'   | .---. .-'  '-.'  .--.; .---. ,--.--.|  .--. '   ",
      "  |  |'.'|  || .-. :'-.  .-'|  |    | .-. ||  .--'|  '--'.'   ",
      "  |  |   |  ||  ---.  |  |  |  '--.;' '-' '|  |   |  |.  .    ",
      "  '--'   '--' `----'  '--'   '-----' '---' '--'   '--' '--'   ")
  ),
  sep = "\n")

  cat(crayon::red(
    c("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
  ),
  sep = "\n")
  cat("\n")

  # Method selection message
  if (identical(method, 1)) {
    cat(crayon::cyan("Method " %+% green$underline$bold ("1") %+% " has been selected\n"))
    cat(crayon::cyan("Used formula: y ~ s(order)\n"))
  } else if (identical(method, 2)) {
    cat(crayon::cyan("Method " %+% green$underline$bold ("2") %+% " has been selected\n"))
    cat(crayon::cyan("Used formula: y ~ s(order, batch)\n"))
  } else {
    cat(crayon::cyan("Method " %+% green$underline$bold ("2") %+% " has been selected by default\n"))
    cat(crayon::cyan("Used formula: y ~ s(order, batch)\n"))
    method <- 2
  }

  # Coerce vectors and data
  order <- as.numeric(order)
  if (method == 2) batch <- as.numeric(batch)

  int_df <- as.data.frame(int_data)
  stopifnot(nrow(int_df) == length(order), nrow(int_df) == length(class))
  if (method == 2) stopifnot(nrow(int_df) == length(batch))

  # log2 transformation
  int_log <- as.data.frame(lapply(int_df, function(x) log2(as.numeric(x) + 1.1)))

  # Fitting & Predictions
  qc_id <- grep(qc_label, class)
  if (length(qc_id) == 0) stop("No QC rows detected with qc_label = ", deparse(qc_label))

  pboptions(style = 1, char = "-")
  message("Fitting GAMs on QC samples...")

  if (method == 1) {
    fit_data_qc <- data.frame(order = order[qc_id], int_log[qc_id, , drop = FALSE])
    fits <- pbapply::pblapply(
      X = seq_len(ncol(int_log)),
      FUN = function(j) mgcv::gam(fit_data_qc[[j + 1]] ~ s(order), data = fit_data_qc, method = "REML")
    )
    message("Predicting for all samples...")
    preds <- pbapply::pblapply(
      X = fits,
      FUN = function(fit) stats::predict(fit, newdata = data.frame(order = order))
    )
  } else {
    fit_data_qc <- data.frame(order = order[qc_id], batch = batch[qc_id], int_log[qc_id, , drop = FALSE])
    fits <- pbapply::pblapply(
      X = seq_len(ncol(int_log)),
      FUN = function(j) mgcv::gam(fit_data_qc[[j + 2]] ~ s(order, batch), data = fit_data_qc, method = "REML")
    )
    message("Predicting for all samples...")
    preds <- pbapply::pblapply(
      X = fits,
      FUN = function(fit) stats::predict(fit, newdata = data.frame(order = order, batch = batch))
    )
  }

  # Correction
  qc_means <- vapply(int_log[qc_id, , drop = FALSE], mean, numeric(1), na.rm = TRUE)
  corrected_log <- mapply(
    FUN = function(col_j, pred_j, qc_mean_j) col_j - pred_j + qc_mean_j,
    col_j = int_log, pred_j = preds, qc_mean_j = qc_means, SIMPLIFY = FALSE
  )
  res <- as.data.frame(lapply(corrected_log, function(x) 2^x))
  rownames(res) <- rownames(int_df)
  colnames(res) <- colnames(int_df)
  res
}
