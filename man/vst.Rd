\name{vst}
\alias{vst}
\title{
Variance Stabilizing Transformation
}

\description{
Perform a Variance Stabilizing Transformation (VST) of a matrix of count data.
}

\usage{
  vst(x, method = "anscombe.nb", lib.size = NULL, cpm = FALSE, dispersion = NULL, design = NULL)
}

\arguments{
  \item{x}{ A matrix of counts. Rows are features (such as genes), and columns are samples.}
  \item{method}{ VST to use, see details. }
  \item{lib.size}{ Optional, estimated if not given. }
  \item{cpm}{ Should the output be in log2 Counts Per Million, rather than simply log2. }
  \item{dispersion}{ Optional, estimated if not given. Dispersion parameter of the negative binomial distribution of the data. }
  \item{design}{ Optional. If dispersion isn't given, a design matrix to use when estimating dispersion. }
}

\details{

Dispersion: 

edgeR's estimate of the common dispersion of the count matrix would be a good choice of dispersion.

If "dispersion" is not given, it is chosen so as to minimize sd(RSS)/mean(RSS) of the Residual Sum of Squares of each feature. Residuals are calculated from the linear model specified by the parameter "design".

If "design" also isn't given, Principal Components Analysis is used and the residuals are calculated by removing all but the smallest axis of variation of the data.

}

\value{
A transformed matrix.
}

\references{

}

\author{
Paul Harrison
}

