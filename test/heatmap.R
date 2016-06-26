
source("test/common.R")

y <- matrix(rnorm(100), nrow=25)
rownames(y) <- as.character(seq_len(nrow(y))**10)
colnames(y) <- as.character(seq_len(ncol(y))**10)

save_ggplot(
    "heatmap",
    varistran::plot_heatmap(y)
)
