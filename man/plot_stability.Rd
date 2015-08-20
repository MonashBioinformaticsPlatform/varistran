\name{plot_stability}
\alias{plot_stability}
\title{
Stability plot.
}
\description{
Produce a ggplot object containing a plot of residual standard deviation against mean count.
}
\usage{
plot_stability(y, x = NULL, design = NULL, bins = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Transformed counts matrix. }
  \item{x}{ Optional, original counts matrix. }
  \item{design}{ Matrix specifying a linear model with which to calculate residuals. }
  \item{bins}{ Number points in the graph. }
}
\details{
Genes are partitioned evenly into "bins" bins by average expression level. Mean residual standard deviation is plotted against mean count.
}
\value{
A ggplot object.

This must be print()-ed to actually plot.
}
\author{
Paul Harrison
}
\examples{
}
