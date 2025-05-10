# 16S-rRNAseq
Analyze 16S microbiome data using dada2

This project performs a comprehensive 16S rRNA microbiome analysis using the DADA2 pipeline in R on a mock community dataset (MiSeq_SOP). The goal is to process raw Illumina paired-end reads to obtain high-resolution amplicon sequence variants (ASVs), assess microbial diversity, and visualize taxonomic composition over time.

Key Steps:

Preprocessing and Quality Control:

Raw forward and reverse reads were loaded and quality profiles were visualized.

Low-quality sequences were filtered and trimmed using filterAndTrim() with parameters to remove reads with ambiguous bases and poor quality tails.

Error Learning and ASV Inference:

Error rates were learned from the filtered reads.

DADA2's core algorithm was used to infer ASVs separately for forward and reverse reads.

Paired reads were merged to form full denoised sequences.

Sequence Table Construction and Chimera Removal:

A sequence table was generated, summarizing ASVs across samples.

Chimeric sequences were identified and removed using a consensus-based approach.

Taxonomic Assignment:

Taxonomy was assigned to each ASV using a pre-trained SILVA database.

Pipeline Tracking:

A tracking table was created to monitor the number of reads retained after each processing step per sample.

Validation Against Mock Community:

The inferred ASVs were compared with known reference sequences from the mock community to assess the accuracy of the pipeline.

Diversity Analysis:

Shannon and Simpson diversity indices were computed to assess alpha diversity across timepoints.

Temporal patterns were visualized using scatter plots colored by sampling period (Early vs Late).

Beta Diversity and NMDS:

Bray-Curtis dissimilarity was calculated and visualized using non-metric multidimensional scaling (NMDS) to explore community structure.

Taxonomic Composition Visualization:

The top 20 abundant ASVs were selected.

A stacked bar plot was used to compare taxonomic distribution between early and late samples, with a custom color scheme based on family-level classification.
