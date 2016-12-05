
# Simulate installation
if (!("installed" %in% commandArgs(trailingOnly=TRUE))) {
    system("development/pre-commit")
    devtools::load_all(export_all=FALSE)
    detach("package:varistran")
}

library("ggplot2")
library("scales")
library("Biobase")
library("biomaRt")

if (!dir.exists("test_output"))
    dir.create("test_output")


load_bottomly <- function() {
    if (!file.exists("test_output/bottomly_eset.RData"))
        download.file(
            "http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData",
            "test_output/bottomly_eset.RData")
    load("test_output/bottomly_eset.RData")

    bottomly.eset
}

load_bottomly_gene_names <- function(bottomly.eset) {
    if (!file.exists("test_output/bottomly_biomart.RData")) {
        #ensembl <- useMart("ensembl")
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
        mouse_ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
        bottomly.bm <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id", values=rownames(bottomly.eset), mart=mouse_ensembl)
        save(bottomly.bm, file="test_output/bottomly_biomart.RData")
    }

    load("test_output/bottomly_biomart.RData")
    bottomly.bm$external_gene_name[
        match(rownames(bottomly.eset),bottomly.bm$ensembl_gene_id)
    ]
}



save_plot <- function(prefix, func, width=4,height=3.5) {
    filename <- paste0("test_output/",prefix,".pdf")
    cat("Plotting",filename,"/ .png\n")

    pdf(filename, width=width, height=height)
    func()
    dev.off()

    filename <- paste0("test_output/",prefix,".png")
    png(filename, width=width, height=height, units="in", res=300)
    func()
    dev.off()
}


save_ggplot <- function(prefix, func, ...) {
    save_plot(prefix, function() print(func), ...)
}


my_theme <-
    theme_minimal() +
    theme(
        panel.border = element_rect(fill=NA),
        panel.grid.minor = element_blank())
