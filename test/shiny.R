

source("test/common.R")

library("Biobase")
library("biomaRt")

bottomly.eset <- load_bottomly()
counts <- exprs(bottomly.eset)
y <- varistran::vst(counts)

print( varistran::shiny_report(y, counts=counts) )
