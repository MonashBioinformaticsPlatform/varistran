# Varistran

Varistran is an R package providing a Variance Stabilizing Transformation appropriate for RNA-Seq data, and a variety of diagnostic plots based on such transformation.

* [Function reference](http://logarithmic.net/varistran/reference/index.html)

* [Online demo](http://rnasystems.erc.monash.edu:3838/pfh/2015/demo-varistran)

* [A slideshow describing Varistran](http://rnasystems.erc.monash.edu:3838/pfh/2016/varistran/)

* [Poster for ABACBS 2015](doc/varistran-poster-abacbs-2015.pdf) [(on F1000, doi: 10.7490/f1000research.1110757.1)](http://f1000research.com/posters/4-1041)

Varistran is developed by Paul Harrison (paul.harrison@monash.edu, [@paulfharrisson](https://twitter.com/paulfharrison)) for the [Monash Bioinformatics platform](https://platforms.monash.edu/bioinformatics/).

## Install

Varistran is most easily installed from GitHub using devtools:

```
install.packages("devtools")
devtools::install_github("MonashBioinformaticsPlatform/varistran")
```

### Dependencies

By default, library size estimation is by TMM, implemented in edgeR from BioConductor. You will need to install this manually if you haven't already:

```
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

Varistran also uses various CRAN packages, which will be installed automatically using the `devtools::install_github` method above. Refer to the DESCRIPTION file for the full list of dependencies.

## Example usage

I intend to keep adding features to Varistran. Therefore, I recommend using its functions via `varistran::` rather than attaching its namespace with `library(varistran)`.

Say you have a count matrix `counts` and a design matrix `design`. To perform a variance stabilizing transformation:

```
y <- varistran::vst(counts, design=design)
```

By default, Anscombe's variance stabilizing transformation for the negative binomial distribution is used. This behaves like log2 for large counts (log2 Counts-Per-Million if `cpm=T` is given).

An appropraite dispersion is estimated with the aid of the design matrix. If omitted, this defaults to a column of ones, for blind estimation of the dispersion. This might slightly over-estimate the dispersion. A third possibility is to estimate the dispersion with edgeR.

### Diagnostic plots

`plot_stability` allows assessment of how well the variance has been stabilized:

```
varistran::plot_stability(y, counts, design=design)
```

`plot_biplot` provides a two-dimensional overview of your samples and genes using Principle Components Analysis (similar concept to `plotMDS` in limma):

```
varistran::plot_biplot(y)
```

<img src="doc/biplot-example.png" height="300">

`plot_heatmap` draws a heatmap.

<img src="doc/heatmap-example.png" height="300">


### Shiny report

Varistran's various diagnostic plots are also available as a Shiny app, which can be launched with:

```
varistran::shiny_report(y, counts)
```

If `y` is not given it is calculated with `varistran::vst(counts)`.

```
varistran::shiny_report(counts=counts)
```

* [Online demo](http://rnasystems.erc.monash.edu:3838/pfh/2015/demo-varistran)


## Test suite

After downloading the source code, a suite of tests can be run with:

```
make test
```

Outputs are places in a directory called `test_output`.


## Links

* [Monash Bioinformatics Platform, Monash University](https://platforms.monash.edu/bioinformatics)

* [RNA Systems Laboratory, Monash University](http://rnasystems.erc.monash.edu)


