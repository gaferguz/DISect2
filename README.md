This package is designed as a wrapper to perform all desired pairwise comparisons using DESeq2 (using Wald test with "ash" shrinkage correction).

Functions:

- Note2GO: turns a protein annotation file (typical eggnog format) into two tables in the format required for GSEA and ORA analysis (term2gene and a term2name).
- DISect2: takes both a sample metadata file and the path to the count files to perform DESeq2 + GSEA (it can handle htseqcounts or kallisto h5), making all possible comparisons for a given vector. A subvector can also be specified (subgroup that is used to split the count matrix and compare on each group of the main vector)

sample table format:
-
.
.
.
