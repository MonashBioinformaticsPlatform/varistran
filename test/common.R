
source("R/varistran.R")
source("R/biplot.R")

library("ggplot2")
library("scales")


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


save_plot <- function(prefix, func, width=4,height=3.5) {
    filename <- paste0("test_output/",prefix,".pdf")
    cat("Plotting",filename,"/ .png\n")

    pdf(filename, width=width, height=height)
    func()
    dev.off()

    filename <- paste0("test_output/",prefix,".png")
    png(filename, width=width, height=height, units="in", res=120)
    func()
    dev.off()

}


save_ggplot <- function(prefix, func, ...) {
    save_plot(prefix, function() print(func), ...)
}


my_theme <- theme_minimal() + theme(panel.border = element_rect(fill=NA), panel.grid.minor = element_blank())