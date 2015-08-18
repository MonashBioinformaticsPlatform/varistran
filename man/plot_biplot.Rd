\name{plot_biplot}
\alias{plot_biplot}
\title{
Biplot of expression data
}
\description{
Produce a ggplot containing a biplot of expression data.
}
\usage{
plot_biplot(x, sample_labels = NULL, balance = 0.25)
}
\arguments{
  \item{x}{ 
      Matrix of expression levels, with features (eg genes) as rows and samples as columns.
  }
  \item{sample_labels}{
      Sample labels.
  }
  \item{balance}{
      Relative scaling of features and samples.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A ggplot object.
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
