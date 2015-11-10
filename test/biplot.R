
source("test/common.R")

bottomly.eset <- load_bottomly()

if (!("gene_names" %in% ls())) {
    ensembl <- useMart("ensembl")
    mouse_ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
    result <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id", values=rownames(bottomly.eset), mart=mouse_ensembl)
    rownames(result) <- result$ensembl_gene_id
    gene_names <- result[rownames(bottomly.eset), "external_gene_name"]
}

counts <- exprs(bottomly.eset)
good <- rowSums(counts) >= 1
counts <- counts[good, ]
strain <- phenoData(bottomly.eset)$strain
experiment.number <- factor( phenoData(bottomly.eset)$experiment.number )
design <- model.matrix(~ strain + experiment.number)

labels <- paste0(substr(strain,1,1), experiment.number)

y <- varistran::vst(counts, design=design)

#y <- varistran::vst(counts, method="anscombe.nb.simple", dispersion=1)

#library(limma)
#y <- voom(counts, lib.size=colSums(counts)*calcNormFactors(counts))$E

save_ggplot(
    "biplot",
    varistran::plot_biplot(
        y,
        sample_labels=labels,
        feature_labels=gene_names[good],
        text_size=0.025
        ) + my_theme,
    width=9,height=9
    )

yy <- y
colnames(yy) <- NULL
rownames(yy) <- NULL
save_ggplot(
    "biplot-unlabelled",
    varistran::plot_biplot(yy, text_size=0.025) + my_theme,
    width=9,height=9
    )
