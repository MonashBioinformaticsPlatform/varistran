#
# Plot several transformations
#

source("test/common.R")

x <- 0:10

dispersion <- 0.1

dfs <- list()

funcs <- list(
    "log(x+c)" = varistran::vst_methods$anscombe.nb.simple$vst,
    "Anscombe" = varistran::vst_methods$anscombe.nb$vst,
    "Naive" = varistran::vst_methods$naive.nb$vst,
    "log(x)" = function(x, dispersion) log2(x)
)

dfs <- list()
for(name in names(funcs)) {
    dfs[[name]] <- data.frame(
        name = rep(name, length(x)),
        x = x,
        y = funcs[[name]](x, dispersion)
        )
}

bigdf <- do.call(rbind, dfs)

bigdf <- bigdf[ is.finite(bigdf$y), ]

bigdf$name <- factor(bigdf$name, names(funcs))

save_ggplot(
  "vsts",
  ggplot(bigdf, aes(x=x,y=y,color=name))
  + geom_line()
  + geom_point(size=2)
  + labs(x="x", y="y", color="Transformation", title="Variance stabilizing transformations\ndispersion=0.1")
  + scale_x_continuous(breaks=(0:5)*2)
  + my_theme
  + theme(legend.justification=c(1,0), legend.position=c(1,0))
)
