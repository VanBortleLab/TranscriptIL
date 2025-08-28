#' Identify nested introns and assign intron clusters
#'
#' @description
#' For each gene in a data frame of introns, classify introns as
#' *Parent*, *Nested*, or *Orphan* based on genomic coordinates.
#' Parent introns contain at least one smaller intron within their span,
#' Nested introns fall completely inside a parent intron, and Orphan introns
#' do not overlap other introns. The function also assigns a cluster number
#' (`Intron_cluster`) to each parent–nested group.
#'
#' @param df A data.frame containing at least the columns:
#'   - `gene`: gene identifier
#'   - `Start`: intron start coordinate (numeric)
#'   - `End`: intron end coordinate (numeric)
#'
#' @return A data.frame with two new columns:
#'   - `Nested`: classification of each intron as "Parent", "Nested", or "Orphan"
#'   - `Intron_cluster`: cluster number for each parent–nested group (orphans remain "X")
#'
#' @details
#' Introns are processed gene by gene. For each intron, overlaps are used
#' to determine whether it contains other introns (Parent), is contained
#' within another (Nested), or stands alone (Orphan). Parent–nested groups
#' are assigned unique cluster numbers for downstream analysis.
#'
#' @examples
#' \dontrun{
#' test_df <- data.frame(
#'   gene  = c("Gene1","Gene1","Gene1","Gene2"),
#'   Start = c(100,150,400,200),
#'   End   = c(500,200,450,300)
#' )
#' nested_intron(test_df)
#' }
#' @export



nested_intron <- function(df){

#create columns to store information
df$Nested = "X"
df$Intron_cluster = "X"
genes = unique(df$gene)
cluster_num = 1
output = NULL

for(z in 1:length(genes)){
  s = subset(df, df$gene == genes[z])

  #### step 1, determine parent and nested introns:
  for(i in 1:nrow(s)){
    tmpP = subset(s, s$Start <= s$Start[i] & s$End >= s$End[i])
    tmpN = subset(s, s$Start >= s$Start[i] & s$End <= s$End[i])

    if(nrow(tmpP) > 1){
      s$Nested[i] = "Nested" ### if there are any parents at all, this intron is "nested"
    }else if(nrow(tmpN) > 1){
      s$Nested[i] = "Parent" ### otherwise, if there are only nested introns, this intron is "parent"
    }else{
      s$Nested[i] = "Orphan" ### otherwise, this intron is an "orphan"
    }
  }

  #### step 2, assign cluster number to each nested intron group:
  parents = subset(s, s$Nested == "Parent")
  nested	= subset(s, s$Nested == "Nested")
  orphans = subset(s, s$Nested == "Orphan")

  if(nrow(parents) > 0){
    for(p in 1:nrow(parents)){
      cluster_num = cluster_num + 1

      parents$Intron_cluster[p] = cluster_num
      tmpN = subset(nested, nested$Start >= parents$Start[p] & nested$End <= parents$End[p])
      tmpN$Intron_cluster = cluster_num
      grouping = rbind(parents[p,], tmpN)
      output = rbind(output, grouping)
    }
    output = rbind(output, orphans)
  }else if(nrow(parents) == 0){
    output = rbind(output, orphans)
  }
}
return(output)
}


