# asinh(ab) / log(2) + const, behaving as log2(a) in the limit as a -> large
.log2ish.asinh <- function(a,b) {
    ab <- a * b
    log2(ab+sqrt(1+ab*ab)) - log2(b) - 1
}

#' @export
vst.methods <- list(
    naive.poisson = list(
        is.logish        = FALSE,
        needs.dispersion = FALSE,
        vst              = function(x) sqrt(x)
    )
    ,
    naive.nb = list(
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            .log2ish.asinh(sqrt(x),sqrt(dispersion)) * 2.0
    )
    ,
    anscombe.poisson = list(
        is.logish        = FALSE,
        needs.dispersion = FALSE,
        vst              = function(x)
            sqrt(x+0.375)
    )
    ,
    anscombe.nb.simple = list(
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            log2(x + 0.5/dispersion)
    )
    ,
    anscombe.nb = list(
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            .log2ish.asinh(sqrt(x+0.375),sqrt(1/(1/dispersion-0.75))) * 2.0
    )
)



.dispersion.score <- function(mat, design=NULL) {
    if (is.null(design))
        design <- matrix(1, ncol=1, nrow=ncol(mat))

    # Residuals, rotated from the subspace orthogonal
    # to the linear model that they inhabit (using QR decomposition)
    residuals <- mat %*% MASS::Null(design)

    # Residual standard deviation
    rsd <- sqrt(rowMeans(residuals*residuals))

    sd(rsd) / mean(rsd)
}


.optimal.dispersion <- function(x, method="anscombe.nb", lib.size=NULL, design=NULL) {
    x <- x[ rowMeans(x) >= 5, ]

    optimize(
         function(d) {
             .dispersion.score(
                 vst(x,method=method,lib.size=lib.size,dispersion=d),
                 design=design)
         },
         lower = 1e-4,
         upper = 1.0
    )$minimum
}




#' Variance Stabilizing Transformation
#'
#' Perform a Variance Stabilizing Transformation (VST) of a matrix of count
#' data.
#'
#'
#' Several methods are available. "anscombe.nb" is recommended.
#'
#' Methods:
#'
#' "anscombe.nb": Default, asinh(sqrt((x+3/8)/(1/dispersion-3/4))). Anscombe's
#' VST for the negative binomial distribution.
#'
#' "anscombe.nb.simple": log(x+0.5/dispersion), a simplified VST also given by
#' Anscombe.
#'
#' "anscombe.poisson": sqrt(x+3/8). Anscombe's VST for the Poisson
#' distribution. Only appropriate if you know there is no biological noise.
#'
#' "naive.nb": asinh(sqrt(x/dispersion)). Resultant variance is slightly
#' inflated at low counts.
#'
#' "naive.poisson": sqrt(x). Resultant variance is slightly inflated at low
#' counts.
#'
#' Dispersion:
#'
#' edgeR's estimate of the common dispersion of the count matrix would be a
#' reasonable choice of dispersion. However Poisson noise in RNA-Seq data may
#' be over-dispersed, in which case a slightly smaller dispersion may work
#' better. I recommend not providing a dispersion and letting varistran pick an
#' appropriate value.
#'
#' If "dispersion" is not given, it is chosen so as to minimize sd(residual
#' s.d.)/mean(residual s.d.). Residuals are calculated from the linear model
#' specified by the parameter "design".
#'
#' If "design" also isn't given, a linear model containing only an intercept
#' term is used. This may lead to an over-estimate of the dispersion, so do
#' give a design if possible.
#'
#' @param x A matrix of counts. Rows are genes (or other features), and columns
#' are samples.
#' @param method VST to use, see details.
#' @param lib.size Optional, estimated if not given.
#' @param cpm Should the output be in log2 Counts Per Million, rather than
#' simply log2.
#' @param dispersion Optional, estimated if not given. Dispersion parameter of
#' the negative binomial distribution of the data.
#' @param design Optional. If dispersion isn't given, a design matrix to use
#' when estimating dispersion.
#' @return A transformed matrix.
#' @author Paul Harrison
#' @references Anscombe, F.J. (1948) "The transformation of Poisson, binomial,
#' and negative-binomial data", Biometrika 35 (3-4): 246-254
#'
#' @export
vst <- function(x, method="anscombe.nb", lib.size=NULL, cpm=FALSE, dispersion=NULL, design=NULL) {
    x <- as.matrix(x)

    method.info <- vst.methods[[method]]
    is.null(method.info) && stop("Unknown method")

    if (is.null(lib.size)) {
        lib.size <- colSums(x) * edgeR::calcNormFactors(x)
    }

    mean.size <- mean(lib.size)

    x.norm <- t(t(x) * (mean.size / lib.size))

    if (method.info$needs.dispersion) {
        if (is.null(dispersion)) {
            dispersion <- .optimal.dispersion(x,method=method,lib.size=lib.size,design=design)
            cat("Dispersion estimated as ",dispersion,"\n",sep="")
        }

        result <- method.info$vst(x.norm, dispersion)
    } else {
        result <- method.info$vst(x.norm)
    }

    if (cpm) {
        method.info$is.logish ||
            stop("Counts Per Million is not meaningful with this transform, use cpm=FALSE")
        result <- result + log2(1e6/mean.size)
    }

    if (!is.null(colnames(x)))
        colnames(result) <- colnames(x)

    if (!is.null(rownames(x)))
        rownames(result) <- rownames(x)

    attr(result, "lib.size") <- lib.size
    attr(result, "cpm") <- cpm
    attr(result, "method") <- method
    if (method.info$needs.dispersion)
        attr(result, "dispersion") <- dispersion

    result
}
