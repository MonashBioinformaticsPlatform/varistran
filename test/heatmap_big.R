
source("test/common.R")

n <- 1000
p <- 5
y <- matrix(rnorm(n*p), nrow=n)

save_ggplot(
    "heatmap_big",
    varistran::plot_heatmap(y),
    width=8,
    height=8)
