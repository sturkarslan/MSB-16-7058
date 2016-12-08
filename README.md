# MSB-16-7058
  Supporting files for manuscript Turkarslan et al., MSB-16-7058

# File Descriptions
## EGRIN-model-cytoscape.cys
  Cytoscape file for visualizing EGRIN Model

## cMonkey-Inferelator.RData
  RData file for cMonkey and Inferelator algorithm runs

## analyzeRNASeq.py
  Custom workflow for analysis of RNAseq data.
  1. Read quality control (FastQC)
  2. Trimming and filtering (Trimmomatic)
  3. Alignment to reference genome (STAR)
  4. Counting reads per transcript (htseq-count)
  5. Calculating FPKMs (Cufflinks)

## rnaseqStats.R
  Custom workflow for analysis of RNAseq data
  1. read stats and quality metrics 
  2. Read count normalizations
  3. Visualizations

## transcriptomeProteomeComparison.py
  Tool to compare correlation between transcriptomics and proteomics data.
  Originally cloned from repository: https://github.com/adelomana/hondatutak

## trendSelector.py
  Tool to find systematically trends of regulation in protein expression data.
  Originally cloned from repository: https://github.com/adelomana/hondatutak
