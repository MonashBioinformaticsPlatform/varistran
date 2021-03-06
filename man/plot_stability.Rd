% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_stability}
\alias{plot_stability}
\title{Stability plot.}
\usage{
plot_stability(y, x = NULL, design = NULL, bins = 20)
}
\arguments{
\item{y}{Transformed counts matrix.}

\item{x}{Optional, original counts matrix.}

\item{design}{Matrix specifying a linear model with which to calculate
residuals.}

\item{bins}{Number points in the graph.}
}
\value{
A ggplot object.

This must be print()-ed to actually plot.
}
\description{
Produce a ggplot object containing a plot of residual standard deviation
against mean count.
}
\details{
Genes are partitioned evenly into "bins" bins by average expression level.
Mean residual standard deviation is plotted against mean count.

If the variance stabilizing transformation has been successful, this plot should be close to a horizontal line. However it is normal for the standard deviation to drop off for counts below 5.
}
\author{
Paul Harrison
}
