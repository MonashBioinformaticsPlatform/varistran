\name{plot_biplot}
\alias{plot_biplot}
\title{
Biplot of expression data
}
\description{
Produce a ggplot object containing a biplot of expression data.
}
\usage{
plot_biplot(x, sample_labels = NULL, balance = 0.25)
}
\arguments{
  \item{x}{ 
      Matrix of expression levels, with features (eg genes) as rows and samples as columns. For example, you could use the output of varistran::vst here.
  }
  \item{sample_labels}{
      Sample labels.
  }
  \item{balance}{
      Relative scaling of features and samples.
  }
}
\details{
Biplot based on the Singular Value Decomposition of the matrix x. The dimensions corresponding to the two largest singular values are shown.

Genes are shown in blue and samples in red.

The dot product of the gene and sample vectors approximates the difference from the average expression level of that gene in that sample.

Sample points (red) are scaled to have the same variance in the two dimensions. Therefore the gene points (blue) may have greater variance along dimension 1 if dimension 1 explains more of the variance than dimension 2.

}
\value{
A ggplot object.

This must be print()-ed to actually plot.
}
\references{
}
\author{
Paul Harrison
}
\note{
}

\seealso{
}
\examples{

# Assuming counts is a matrix of read counts.
y <- varistran::vst(counts)
print( varistran::plot_biplot(y) )

}
