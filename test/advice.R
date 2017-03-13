source("test/common.R")

bottomly.eset <- load_bottomly()

gene_names <- load_bottomly_gene_names(bottomly.eset)

counts <- exprs(bottomly.eset)
strain <- phenoData(bottomly.eset)$strain
experiment.number <- factor( phenoData(bottomly.eset)$experiment.number )
design <- model.matrix(~ strain + experiment.number)

y <- varistran::vst(counts, design=design)

print(varistran::vst_advice(y))
