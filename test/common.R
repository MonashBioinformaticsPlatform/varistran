
source("R/varistran.R")

library("ggplot2")


if (!dir.exists("test_output"))
    dir.create("test_output")


save.plot <- function(prefix, func, width=4,height=3.5) {
    filename <- paste0("test_output/",prefix,".pdf")
    cat("Plotting",filename,"\n")
    pdf(filename, width=width, height=height)
    func()
    dev.off()
}


save.ggplot <- function(prefix, func, ...) {
    save.plot(prefix, function() print(func), ...)
}


my_theme <- theme_minimal() + theme(panel.border = element_rect(fill=NA), panel.grid.minor = element_blank())