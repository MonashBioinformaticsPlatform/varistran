# Varistran

Varistran is an R package providing a Variance Stabilizing Transformation appropriate for RNA-Seq data, and a variety of diagnostic plots based on such transformation.

* [Online demo](http://rnasystems.erc.monash.edu:3838/pfh/2015/demo-varistran)

Varistran is developed by Paul Harrison (paul.harrison@monash.edu, [@paulfharrisson](https://twitter.com/paulfharrison)) for the [Monash Bioinformatics platform](https://platforms.monash.edu/bioinformatics/).

* [Documentation (PDF)](http://rnasystems.erc.monash.edu/doc/varistran.pdf)

* [Poster for ABACBS 2015](doc/varistran-poster-abacbs-2015.pdf) [(on F1000, doi: 10.7490/f1000research.1110757.1)](http://f1000research.com/posters/4-1041)

## Install

Varistran is most easily installed from GitHub using devtools:

```
install.packages("devtools")

devtools::install_github("MonashBioinformaticsPlatform/varistran")
```

### Dependencies

Various diagnostic plotting functions in Varistran use `grid`, `ggplot2` and `ggdendro`. Varistran can also produce an interactive `shiny` report. `seriation` is used for heatmaps.

```
install.packages(c("grid","ggplot2","ggdendro","seriation","shiny"))
```

By default, library size estimation is by TMM, implemented in edgeR from BioConductor:

```
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

## Example usage

I expect to keep adding features to Varistran. Therefore, I recommend using its functions via `varistran::` rather than attaching its namespace with `library(varistran)`.

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

![Example of a biplot](doc/biplot-example.png)

### Shiny report

Varistran's various diagnostic plots are also available as a Shiny app, which can be launched with:

```
varistran::shiny_report(y, counts)
```

* [Online demo](http://rnasystems.erc.monash.edu:3838/pfh/2015/demo-varistran)


## Links

* [Monash Bioinformatics Platform, Monash University](https://platforms.monash.edu/bioinformatics)

* [RNA Systems Laboratory, Monash University](http://rnasystems.erc.monash.edu)
