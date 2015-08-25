# Varistran

An R package providing Variance Stabilizing Transformations appropriate for RNA-Seq data.

* [Documentation (PDF)](http://rnasystems.erc.monash.edu/doc/varistran.pdf)

## Install

Varistran is most easily installed from GitHub using devtools:

```
install.packages("devtools")

devtools::install_github("MonashBioinformaticsPlatform/varistran")
```

### Dependencies

Various diagnostic plotting functions in Varistran use `ggplot2`:

```
install.packages("ggplot2")
```

By default, library size estimation is by TMM, implemented in edgeR from BioConductor:

```
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

## Example usage

Say you have a count matrix `counts` and a design matrix `design`. To perform a variance stabilizing transformation:

```
y <- varistran::vst(counts, design=design)
```

By default, Anscombe's variance stabilizing transformation for the negative binomial distribution is used. This behaves like log2 for large counts (log2 Counts-Per-Million if `cpm=T` is given).

An appropraite dispersion is estimated with the aid of the design matrix. (If omitted, this defaults to a column of ones, for blind estimation of the dispersion. This might slightly over-estimate the dispersion. A third possibility is to estimate the dispersion with edgeR.)

### Diagnostic plots

`plot_stability` allows assessment of how well the variance has been stabilized:

```
varistran::plot_stability(y, counts, design=design)
```

`plot_biplot` provides a two-dimensional overview of your samples and genes using Principle Components Analysis (similar to `plotMDS` in limma):

```
varistran::plot_biplot(y)
```
