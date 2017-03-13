---
title: "Varistran: Anscombe's variance stabilizing transformation for RNA-seq gene expression data"
tags:
  - RNA-Seq
  - gene expression
  - variance stabilizing transformation
  - R package
authors:
  - name: Paul Francis Harrison
    orcid: 0000-0002-3980-268X
    affiliation: 1
affiliations:
  - name: Monash Bioinformatics Platform, Monash University
    index: 1
date: 8 March 2017
bibliography: paper.bib
---

# Summary

RNA-seq measures RNA expression levels in a biological sample using high-throughput cDNA sequencing, producing counts of the number of reads aligning to each gene. Noise in RNA-seq read count data is commonly modelled as following a negative binomial distribution, where the variance is a quadratic function of the expression level. However many statistical, machine learning, and visualization methods work best when the noise in a data set has equal variance. Varistran is an R package that uses Anscombe's [-@Anscombe1948] variance stabilizing transformation for the negative binomial distribution to transform RNA-seq count data, so that the noise has equal variance across all measured gene expression levels. The transformed data may be treated as log~2~ transformed gene expression levels, but with variability reduced at low read counts. Varistran also includes a function to open a Shiny report with simple diagnostic visualizations, including a plot to assess how effective the variance stabilization was, a biplot of samples and genes, and a heatmap. This allows defective samples, sample mislabling, and batch effects to be easily identified.

# References

