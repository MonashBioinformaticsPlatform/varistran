
source("test/common.R")

set.seed(563)
n <- 5000
p <- 10
y <- matrix(rnorm(n*p), nrow=n)

save_ggplot(
    "heatmap_big",
    varistran::plot_heatmap(y, show_tree=FALSE),
    width=8,
    height=8)
