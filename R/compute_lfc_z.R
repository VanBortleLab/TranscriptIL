#' Compute log2 fold change (experiment over control) and transform into Z-score
#'
#' @description
#' Given a merged IRFinder table (with replicate means already computed),
#' calculate log2 fold change for `{experiment_name}_{metric}_avg` over
#' `{control_name}_{metric}_avg`, and add a standardized Z-score of that LFC.
#'
#' @param df Data frame (typically output of `compute_group_means()`).
#' @param experiment_name Character prefix of the experiment sample
#'   (e.g., "HEK293T_snaR").
#' @param control_name Character prefix of the control sample
#'   (e.g., "HEK293T_scramble").
#' @param new_name Character prefix of the new name
#'   (e.g., "snaR").
#' @param metric Which averaged metric to use; default `"IRratio"`.
#'   The function expects columns named `{sample_name}_{metric}_avg`.
#' @param pseudocount Small numeric added to numerator and denominator to avoid
#'   log of zero (default `1e-6`).
#'
#' @return The input data frame with two new columns:
#' - `{new_name}_{metric}_lfc`
#' - `{new_name}_{metric}_z`
#'
#' @examples
#' \dontrun{
#' df2 <- compute_lfc_z(
#'   df,
#'   experiment_name = "HEK293T_snaR",
#'   control_name    = "HEK293T_scramble",
#'   new_name = "snaR",
#'   metric = "IRratio"
#' )
#' }
#' @export
compute_lfc_z <- function(df,
                          control_name,
                          experiment_name,
                          new_name,
                          metric = "IRratio",
                          pseudocount = 1e-6) {
  exp_avg_col <- paste0(experiment_name, "_", metric, "_avg")
  ctl_avg_col <- paste0(control_name, "_", metric, "_avg")

  if (!exp_avg_col %in% names(df))
    stop("Missing experiment mean column: ", exp_avg_col)
  if (!ctl_avg_col %in% names(df))
    stop("Missing control mean column: ", ctl_avg_col)

  lfc_col <- paste0(new_name,"_", metric, "_lfc")
  z_col   <- paste0(new_name, "_", metric, "_z")

  lfc <- log2((df[[exp_avg_col]] + pseudocount) / (df[[ctl_avg_col]] + pseudocount))
  z   <- as.numeric(scale(lfc))  # center/scale to Z across rows

  df[[lfc_col]] <- lfc
  df[[z_col]]   <- z
  df
}

