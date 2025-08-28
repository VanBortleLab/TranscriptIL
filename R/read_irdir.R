#' Read IRFinder replicate files from a directory
#'
#' @description
#' Reads one or more IRFinder replicate result files from a directory,
#' constructs a unique intron index (`Chr_Start_End_Name`), and returns
#' a data frame with per-replicate IRratio and SpliceExact values.
#'
#' @param dir Path to a directory containing IRFinder result files
#'   (e.g., IRFinder-IR-dir-1.txt, IRFinder-IR-dir-2.txt).
#' @param sample_name Character prefix to use for naming replicate columns.
#' @param meta_cols Character vector of required metadata columns
#'   (default: c("Chr","Start","End","Name")).
#' @param keep_extra_meta Logical; if TRUE, keep any extra metadata
#'   columns from the first file. Default = FALSE.
#'
#' @return A data.frame with:
#' - `index` column (unique intron identifier)
#' - IRratio and SpliceExact columns for each replicate, named
#'   `{sample_name}_IRratio_repX`, `{sample_name}_SpliceExact_repX`.

#'
#' @examples
#' \dontrun{
#' df <- read_irfinder("/path/HEK293T_SNAR", sample_name = "HEK293T_snaR")
#' head(df)
#' }
#' @export



read_irfinder <- function(dir,
                          sample_name,
                          meta_cols = c("Chr","Start","End","Name"),
                          keep_extra_meta = FALSE) {
  files <- list.files(dir,full.names = TRUE)

  stopifnot(is.character(sample_name), length(sample_name) == 1)

  # choose fast & robust reader if available
  fread <- getNamespace("data.table")$fread
  have_dt <- !is.null(fread)

  read_one <- function(fp) {
    if (have_dt) {
      df <- data.table::fread(fp, sep = "\t", header = TRUE, data.table = FALSE)
    } else {
      df <- utils::read.table(fp, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
    }
    missing_meta <- setdiff(meta_cols, names(df))
    if (length(missing_meta)) {
      stop(sprintf("File '%s' is missing required columns: %s",
                   fp, paste(missing_meta, collapse = ", ")))
    }
    needed <- c(meta_cols, "IRratio", "SpliceExact")
    missing_needed <- setdiff(needed, names(df))
    if (length(missing_needed)) {
      stop(sprintf("File '%s' is missing IRFinder columns: %s",
                   fp, paste(missing_needed, collapse = ", ")))
    }

    # build index
    df$index <- paste(df[[meta_cols[1]]], df[[meta_cols[2]]], df[[meta_cols[3]]], df[[meta_cols[4]]], sep = "_")

    # keep only index + two metrics (+ optional extras from first file)
    keep <- c("index", "IRratio", "SpliceExact")
    if (keep_extra_meta) {
      # keep any additional non-duplicative columns (first file only later)
      extra <- setdiff(names(df), keep)
      keep <- c(keep, extra)
    }
    df[, keep, drop = FALSE]
  }

  lst <- lapply(files, read_one)

  # Merge by index; take meta from the first file
  base <- lst[[1]]
  base_metrics <- setdiff(names(base), c(meta_cols, "index"))
  names(base)[names(base) %in% c("IRratio","SpliceExact")] <-
    paste0(sample_name, c("_IRratio_rep1","_SpliceExact_rep1"))

  if (length(lst) >= 2) {
    for (i in 2:length(lst)) {
      df <- lst[[i]][, c("index","IRratio","SpliceExact"), drop = FALSE]
      names(df)[names(df) == "IRratio"]     <- paste0(sample_name, "_IRratio_rep", i)
      names(df)[names(df) == "SpliceExact"] <- paste0(sample_name, "_SpliceExact_rep", i)
      base <- merge(base, df, by = "index", all = TRUE, sort = FALSE)
    }
  }

  # ensure meta columns are first
  ord <- c( "index",
           grep(sprintf("^%s_", sample_name), names(base), value = TRUE))
  base <- base[, ord, drop = FALSE]



  # sanity checks
  if (anyNA(base$index)) stop("Index contains NA; check input files.")
  if (anyDuplicated(base$index)) warning("Duplicate indices found after merge.")

  return(base)
}


#' Read and merge IRFinder outputs for control and experiment
#'
#' @description
#' Reads replicate IRFinder outputs from two directories (control and experiment),
#' merges by intron index, adds a `gene` column, filters introns based on
#' minimum SpliceExact across all replicates, and drops SpliceExact columns.
#'
#' @param control_dir Path to directory with IRFinder outputs for control sample.
#' @param experiment_dir Path to directory with IRFinder outputs for experiment sample.
#' @param control_name Character prefix for naming control replicate columns.
#' @param experiment_name Character prefix for naming experiment replicate columns.
#' @param splice_min Numeric threshold (default = 10). Introns with minimum
#'   SpliceExact < threshold across all replicates are removed.
#'
#' @return A data.frame containing:
#' - `gene` as the first column
#' - `index` and coordinate metadata
#' - IRratio replicate columns for both control and experiment
#' - Introns failing the SpliceExact filter removed
#' - with chromosome, intron start, intron end information
#'
#' @examples
#' \dontrun{
#' df <- read_irdir(
#'   control_dir     = "/path/HEK293T_SCRAMBLE",
#'   experiment_dir  = "/path/HEK293T_SNAR",
#'   control_name    = "HEK293T_scramble",
#'   experiment_name = "HEK293T_snaR",
#'   splice_min      = 10
#' )
#' head(df)
#' }
#' @export




read_irdir <- function(control_dir,
                            experiment_dir,
                            control_name,
                            experiment_name,
                            splice_min =10) {
  # 1) read replicate files from each directory
  ctrl <- read_irfinder(control_dir,    sample_name = control_name)
  exp  <- read_irfinder(experiment_dir, sample_name = experiment_name)

  # 2) join by intron index
  merged <- merge(ctrl, exp, by="index")

  # 3) add gene column
  fourth <- vapply(strsplit(merged$index, "_"), `[`, character(1), 4)
  genes <- vapply(strsplit(fourth, "/"), `[`, character(1), 1)
  merged$gene <- genes

  # 4) compute row-wise minimum across ALL SpliceExact columns
  splice_cols <- grep("SpliceExact", names(merged), value = TRUE)
  if (!length(splice_cols)) stop("No SpliceExact columns found in merged table.")
  # row-wise min, ignoring NAs
  merged$SpliceExact_min <- apply(merged[, splice_cols, drop = FALSE], 1, function(x) {
    if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  })

  # 5) filter: keep rows with min SpliceExact >= threshold
  keep <- !is.na(merged$SpliceExact_min) & merged$SpliceExact_min >= splice_min
  merged <- merged[keep, , drop = FALSE]

  # 6) drop all SpliceExact columns (including the *_min helper)
  drop_cols <- grep("SpliceExact", names(merged), value = TRUE)
  merged <- merged[, !(names(merged) %in% drop_cols), drop = FALSE]



  #7) add chromosome, intron start, and intron end
  positions = matrix(unlist(strsplit(merged$index, "_")), byrow=T, ncol=4)
  positions = data.frame(positions)

  merged$Chr = positions$X1
  merged$Start = as.integer(positions$X2)
  merged$End = as.integer(positions$X3)

  # 8) move gene to the first column
  front_cols <- c("gene", "Chr", "Start", "End")
  rest_cols  <- setdiff(names(merged), front_cols)
  merged     <- merged[, c(front_cols, rest_cols), drop = FALSE]

  return(merged)
}



