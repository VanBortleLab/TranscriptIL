#' Compute replicate means for control and experiment groups
#'
#' @description
#' Given a merged IRFinder data frame or a user defined data frame, compute mean IRratio across replicates
#' for both control and experiment groups. The function looks for all replicate
#' columns matching `{sample_name}_IRratio_repX`.
#'
#' @param df Data frame produced by `read_irdir()` or user provied data frame(with replicate IRratio columns).
#' @param control_name Prefix for control sample (e.g., "HEK293T_scramble").
#' @param experiment_name Prefix for experiment sample (e.g., "HEK293T_snaR").
#' @param metric Which metric to average (default = "IRratio").
#'
#' @return A data frame with two new columns:
#' - `{control_name}_{metric}_avg`
#' - `{experiment_name}_{metric}_avg`
#' - remove introns with IR ratio of minimum = 0 or maximum = 1  in either group
#'
#' @examples
#' \dontrun{
#' df <- compute_group_means(df, control_name = "HEK293T_scramble",
#'                                experiment_name = "HEK293T_snaR")
#' head(df)
#' }
#' @export
compute_group_means <- function(df,
                                control_name,
                                experiment_name,
                                metric = "IRratio") {
  # pattern to find replicate columns
  ctrl_pattern <- sprintf("^%s_%s_rep\\d+$", control_name, metric)
  exp_pattern  <- sprintf("^%s_%s_rep\\d+$", experiment_name, metric)

  # identify columns
  ctrl_cols <- grep(ctrl_pattern, names(df), value = TRUE)
  exp_cols  <- grep(exp_pattern,  names(df), value = TRUE)

  if (length(ctrl_cols) == 0)
    stop("No replicate columns found for control: ", control_name)
  if (length(exp_cols) == 0)
    stop("No replicate columns found for experiment: ", experiment_name)

  # compute row means
  df[[paste0(control_name, "_", metric, "_avg")]] <-
    rowMeans(df[, ctrl_cols, drop = FALSE], na.rm = TRUE)

  df[[paste0(experiment_name, "_", metric, "_avg")]] <-
    rowMeans(df[, exp_cols, drop = FALSE], na.rm = TRUE)

  # compute min and max per group
  min_ctrl <- apply(df[, ctrl_cols, drop = FALSE], 1, min, na.rm = TRUE)
  max_ctrl <- apply(df[, ctrl_cols, drop = FALSE], 1, max, na.rm = TRUE)

  min_exp  <- apply(df[, exp_cols, drop = FALSE], 1, min, na.rm = TRUE)
  max_exp  <- apply(df[, exp_cols, drop = FALSE], 1, max, na.rm = TRUE)

  # filter: remove introns with 0 min or 1 max in either group
  keep <- min_ctrl != 0 & max_ctrl != 1 &
    min_exp  != 0 & max_exp  != 1

  df <- df[keep, , drop = FALSE]

  df
}


