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