
source("test/common.R")

counts <- exprs(bottomly.eset)
counts <- counts[rowSums(counts) >= 1, ]
strain <- phenoData(bottomly.eset)$strain
experiment.number <- factor( phenoData(bottomly.eset)$experiment.number )
design <- model.matrix(~ strain + experiment.number)

labels <- paste0(substr(strain,1,1), experiment.number)

y <- vst(counts, design=design)

#y <- vst(counts, method="anscombe.nb.simple", dispersion=1)

#library(limma)
#y <- voom(counts, lib.size=colSums(counts)*calcNormFactors(counts))$E

save_ggplot(
    "biplot",
    plot_biplot(y, sample_labels=labels) + my_theme,
    width=6,height=6
    )

