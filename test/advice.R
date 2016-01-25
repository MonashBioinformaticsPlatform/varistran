
source("test/common.R")

bottomly.eset <- load_bottomly()

gene_names <- load_bottomly_gene_names(bottomly.eset)

counts <- exprs(bottomly.eset)

y <- varistran::vst(counts, design=design)

print(varistran::vst_advice(y))