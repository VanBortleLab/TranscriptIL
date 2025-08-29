
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
