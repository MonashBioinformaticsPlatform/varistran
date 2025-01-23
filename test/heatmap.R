
source("test/common.R")

y <- matrix(rnorm(100, mean=5), nrow=25)
rownames(y) <- as.character(seq_len(nrow(y))**10)
colnames(y) <- as.character(seq_len(ncol(y))**10)

save_ggplot(
    "heatmap",
    varistran::plot_heatmap(y)
)

save_ggplot(
    "heatmap_variant_1",
    varistran::plot_heatmap(y, baseline_to=0, scale_to=3, show_tree=FALSE)
)

save_ggplot(
    "heatmap_variant_2",
    varistran::plot_heatmap(y, show_tree=FALSE, show_baseline=FALSE)
)
