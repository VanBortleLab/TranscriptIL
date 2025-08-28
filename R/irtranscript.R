#' Transcript-level intron retention scoring by permutation
#'
#' @description
#' For each gene, this function aggregates intron-level scores (e.g., Z-scores of log2FC(IRratio) or raw IRratio)
#' into a transcript-level "distribution" score. It accounts for nested introns by taking
#' the median of each parent–nested cluster rather than treating every intron independently.
#' The observed median score is compared to a null distribution of medians obtained from
#' random intron sampling, allowing estimation of empirical p-values and Z-scores.
#'
#' @param df A data.frame of introns, must include:
#'   - `gene`: gene identifier
#'   - `Nested`: intron classification ("Parent","Nested","Orphan")
#'   - `Intron_cluster`: cluster ID for parent–nested groups
#'   - at least one column of numeric scores (e.g., Z-scores of log2FC(IRratio) or raw IRratio)
#' @param score Character vector of column names in `df` to analyze (e.g., `"SNAR_IRratio_z"`).
#'   Each score column is processed in turn.
#'@param save_dir User defined directory for saving the final output; one TSV per score will be saved.
#' @param nperms Number of permutations to build the null distribution (default = 100,000).
#'
#' @return A data.frame with one row per gene, including:
#' - `gene`: gene ID
#' - `obs.IRratio`: observed transcript-level median score
#' - `exp.IRratio`: mean expected score from permutations
#' - `pval_low`: empirical p-value for observed ≤ random
#' - `pval_high`: empirical p-value for observed ≥ random
#' - `sd`: standard deviation of permutation distribution
#' - `z_score`: standardized Z of observed vs expected
#' - `num_intron`: number of introns considered for the gene
#' - `num_intron_without_nested`: number of orphans + clusters used
#'
#' @details
#' The algorithm:
#' 1. For each gene, classify introns as orphan or part of clusters.
#' 2. Compute the median score per cluster and per orphan intron.
#' 3. Combine these into an observed transcript-level score.
#' 4. Compare to a null distribution of random medians sampled from all introns.
#' 5. Report p-values, expected score, and Z-score for each gene.


#' @examples
#' \dontrun{
#' # assuming df already has intron-level z-scores
#' res <- irtranscript(df, score = "SNAR_IRratio_z")
#' head(res)
#' }
#' @export

irtranscript <- function(df,score,save_dir,nperms = 100000){
for(zs in score){
  print(zs)
  score = df[,colnames(df)==zs]
  df$score = as.numeric(score)
  f = sort(unique(df$gene))

  if (!all(c("Nested", "Intron_cluster") %in% names(df))) {
    stop("Input df must contain 'Nested' and 'Intron_cluster' columns. ",
         "Run nested_intron(df) first.")
  }

  output = NULL

  for(j in 1:length(f)){
    fs = subset(df, df$gene == f[j])

    ### if there are nested introns, take that cluster's median score:
    intron_orphans = subset(fs, fs$Nested == "Orphan")
    intron_clusters = subset(fs, fs$Nested != "Orphan")

    clusters = unique(intron_clusters$Intron_cluster)
    cluster_med = vector()
    for(c in 1:length(clusters)){
      tmp = subset(fs, fs$Intron_cluster == clusters[c])
      cluster_med = c(cluster_med, stats::median(tmp$score))
    }

    if(length(clusters) == 0){
      obs = stats::median(fs$score)
      num = nrow(fs)
    }else{
      obs = stats::median(c(intron_orphans$score, cluster_med))
      num = length(c(intron_orphans$score, cluster_med))
    }


    nperms = nperms
    perms = NULL
    for(i in 1:nperms){
      rand = sample(df$score, num)
      perms = rbind(perms, stats::median(rand))
    }

    exp = mean(perms)
    pval_l = length(subset(perms, perms <= obs)) / nperms
    pval_h = length(subset(perms, perms >= obs)) / nperms
    sdv = stats::sd(perms)
    z   = (obs - exp) / sdv

    info = cbind(f[j], obs, exp, pval_l, pval_h, sdv, z, nrow(fs), num)

    if(pval_l != 0){
      output = rbind(output, info)
    }else if(pval_l == 0){
      nperms = nperms
      perms = NULL
      for(i in 1:nperms){
        rand = sample(df$score, nrow(fs))
        perms = rbind(perms, stats::median(rand))
      }

      exp = mean(perms)
      pval_l = length(subset(perms, perms <= obs)) / nperms
      pval_h = length(subset(perms, perms >= obs)) / nperms
      sdv = stats::sd(perms)
      z   = (obs - exp) / sdv

      info = cbind(f[j], obs, exp, pval_l, pval_h, sdv, z, nrow(fs), num)
      output = rbind(output, info)
    }
  }

  output = data.frame(output)

  for(i in 2:ncol(output)){
    output[,i] = as.numeric(output[,i])
  }

 colnames(output) <- c("gene","obs.IRratio","exp.IRratio","pval_low","pval_high","sd","z_score","num_intron","num_intron_without_nested")
 utils::write.table(output, file=paste0(save_dir, "/",zs, "_irTranscript.txt"), sep="\t", quote=F, row.names=F)
}

}



