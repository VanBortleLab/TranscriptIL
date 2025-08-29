
<p align="center">
  <img width="752" height="358" alt="Screenshot 2025-08-29 at 11 49 05 AM"
       src="https://github.com/user-attachments/assets/d94f2ace-ab4a-4d5c-91e0-dab35c45ce0b" />
</p>


# Overview
**TranscriptIL** provides transcript-level analysis of intron retention (IR) with replicate-aware summaries, nested intron clustering, and per-gene permutation testing to quantify transcript-level IR changes with associated p-values. The framework also supports intron-level comparisons between control and treatment conditions.

# Features
- **Intron classification (Parent / Nested / Orphan)**: Intron annotations are often complex, with many introns nested within larger “parent” introns. Such structures can significantly impact the statistics of transcript-level analyses. TranscriptIL addresses this by classifying introns based on their genomic coordinates into Parent, Nested, or Orphan categories, ensuring these relationships are properly accounted for in downstream analysis.
- **Transcript-level IR statistics via permutation tests (empirical p, z)**: TranscriptIL quantifies intron retention at the transcript level by collapsing intron-level signals into a gene-level statistic. A permutation-based framework is then applied to generate empirical p-values and z-scores, providing a rigorous assessment of the significance of intron retention changes.
- Simple, pipeable functions to build a full Transcript-level IR analysis in a few lines

# Installation
```r
#Install TranscriptIL from GitHub
remotes::install_github("VanBortleLab/TranscriptIL")
```

# Example
```r
##Core Workflow Example

# Demo input dataframe
df <- data.frame(
  gene  = c("G1","G1","G2"),
  Start = c(100,150,200),
  End   = c(200,180,300),
  index = c("chr1_100_200_G1/ENSG1",
            "chr1_150_180_G1/ENSG1",
            "chr1_200_300_G2/ENSG2"),
  Condition1_IRratio = c(0.3, 0.2, 0.5),
  Condition2_IRratio = c(0.4, 0.3, 0.6)
)

# Step 1: Classify introns (Parent / Nested / Orphan)
df_nested <- nested_intron(df)

# Step 2: Run transcript-level permutation test
tmpdir <- tempdir()
res <- irtranscript(df_nested, 
                    score    = c("Condition1_IRratio","Condition2_IRratio"), 
                    save_dir = tmpdir,
                    nperms   = 1000)

# View result
print(res)
```
# Functions Overview
## core functions
-- `nested_intron(df)`
- Description: Classifies introns as Parent, Nested, or Orphan based on coordinates; assigns intron clusters.
- Parameters:
`df`: Dataframe with `gene`, `Start`, and `End`.
- Returns: Original dataframe with new `Nested` and `Intron_cluster` columns.

-- `irtranscript(df, score, save_dir, nperms = 100000,, seed = Null)`
- Description: Transcript-level analysis via permutation test; treat nested introns as a cluster to perform the permutation.
- Parameters:<br>
`df`: Dataframe (must include Nested and Intron_cluster). <br>
`score`: Column(s) contains IR ratios to test. <br>
`save_dir`: Directory where results are written. <br>
`nperms`: Number of permutations (default = 100000). <br>
`seed`: Random seed (optional, for reproducibility).<br>
- Returns: Dataframe with per-gene statistics: <br>
`gene`: gene information. <br>
`obs.IRratio`: observed transcript-level median score. <br>
`exp.IRratio`: mean expected score from permutations. <br>
`pval_low`: empirical p-value for observed ≤ random. <br>
`pval_high`: empirical p-value for observed ≥ random. <br>
`sd`: standard deviation of observed median score from permutation distribution. <br>
`z_score`: standardized Z of observed vs expected. <br>
`num_intron`: number of introns before considering nested introns. <br>
`num_intron_without_nested`: number of introns without nested introns. <br>


## additional functions
These functions support reading `IRratio` values from `IRFinder` outputs and is designed to calculate the log2 fold change of experimental IRratio over control, providing a normalized measure of intron retention differences between conditions.

-- `read_irdir(control_dir, experiment_dir, control_name, experiment_name, splice_min = 10)` 
- Description: Reads results from all replicates for control and experiment groups produced by IRfinder(e.g. IRFinder-IR-dir-1.txt, IRFinder-IR-dir-2.txt, IRFinder-IR-dir-3.txt ), merges by intron index, and filters low-count introns. <br>
- Parameters: <br>
`control_dir`: Path to control IRFinder files. <br>
`experiment_dir`: Path to experiment IRFinder files. <br>
`control_name`: Label for the control group. <br>
`experiment_name`: Label for the experiment group. <br>
`splice_min`: Minimum SpliceExact threshold (default = 10). <br>
- Returns: Dataframe with `gene`, `index`, `Chr`, `Start`, `End`, and replicate IRratio values for both groups.<br>

-- `compute_group_means(df, control_name, experiment_name, metric = "IRratio")`
- Description: Calculates replicate mean values for control and experiment. <br>
- Parameters: <br>
`df`: Dataframe from read_irdir. <br>
`control_name`: Name of control group (used in column prefix). <br>
`experiment_name`: Name of experiment group (used in column prefix). <br>
`metric`: Column to average (default = "IRratio"). <br>
- Returns: Dataframe with new `{group}_{metric}_avg` columns. <br>







