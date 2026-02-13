This package is thought as a wrapper to perform all desired pairwise comparison in DESeq2 (Wald test implementation with "ash" shrinkage).

Main functions incorporated

- Note2GO: turns protein annotation file (typical eggnog format) to term2gene and term2name needed for GSEA and ORA analysis
- DISect2: takes sample table file and counts directory and performe DESeq2 + GSEA (htseqcounts or h5 kallisto files) and performes all possible comparision for a vector of factors or subfactors (subgroup within the main factor that is used to split de count matrix)

sample table format
