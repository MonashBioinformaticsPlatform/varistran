

source("test/common.R")

library("Biobase")
library("biomaRt")

bottomly.eset <- load_bottomly()
bottomly.gene.names <- load_bottomly_gene_names(bottomly.eset)
counts <- exprs(bottomly.eset)

strain <- phenoData(bottomly.eset)$strain
experiment.number <- factor( phenoData(bottomly.eset)$experiment.number )
design <- model.matrix(~ strain + experiment.number)
labels <- paste0(substr(strain,1,1), experiment.number)


print( 
    varistran::shiny_report(
        #y=varistran::vst(counts, method="anscombe.poisson"),
        #y=varistran::vst(counts)[,],
        counts=counts, 
        sample_labels=labels, 
        feature_labels=bottomly.gene.names))
