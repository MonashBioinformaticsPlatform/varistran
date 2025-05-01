
# Imports
#
#' @importFrom grDevices axisTicks dev.off pdf png postscript
#' @importFrom stats as.dendrogram dist hclust is.leaf optimize order.dendrogram sd
#' @importFrom graphics par plot.new
#' @importFrom Matrix colSums t

# asinh(ab) / log(2) + const, behaving as log2(a) in the limit as a -> large
.log2ish_asinh <- function(a,b) {
    ab <- a * b
    log2(ab+sqrt(1+ab*ab)) - log2(b) - 1
}

vst_methods <- list(
    naive.poisson = list(
        description      = "Naive Poisson Variance Stabilizing Transformation.",
        units            = "sqrt",
        is.logish        = FALSE,
        needs.dispersion = FALSE,
        vst              = function(x) sqrt(x)
    )
    ,
    naive.nb = list(
        description      = "Naive negative binomial Variance Stabilizing Transformation.",
        units            = "log2",
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            .log2ish_asinh(sqrt(x),sqrt(dispersion)) * 2.0
    )
    ,
    anscombe.poisson = list(
        description      = "Anscombe's Variance Stabilizing Transformation for the Poisson distribution.",
        units            = "sqrt",
        is.logish        = FALSE,
        needs.dispersion = FALSE,
        vst              = function(x)
            sqrt(x+0.375)
    )
    ,
    anscombe.nb.simple = list(
        description      = "Anscombe's simplified Variance Stabilizing Transformation for the negative binomial distribution.",
        units            = "log2",
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            log2(x + 0.5/dispersion)
    )
    ,
    anscombe.nb = list(
        description      = "Anscombe's Variance Stabilizing Transformation for the negative binomial distribution.",
        units            = "log2",
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            .log2ish_asinh(sqrt(x+0.375),sqrt(1/(1/dispersion-0.75))) * 2.0
    )
)



.dispersion_score <- function(mat, design=NULL) {
    if (is.null(design))
        design <- matrix(1, ncol=1, nrow=ncol(mat))

    # Residuals, rotated from the subspace orthogonal
    # to the linear model that they inhabit (using QR decomposition)
    residuals <- mat %*% MASS::Null(design)

    # Residual standard deviation
    rsd <- sqrt(rowMeans(residuals*residuals))

    sd(rsd) / mean(rsd)
}


.optimal_dispersion <- function(x, method="anscombe.nb", lib.size=NULL, design=NULL) {
    x <- x[rowMeans(x) >= 5,,drop=FALSE]

    if (nrow(x) == 0 || ncol(x) == 0) {
        warning("Insufficient data to estimate dispersion, default value used.")
        return(1.0)
    }

    optimize(
         function(d) {
             .dispersion_score(
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
#' @examples
#'
#' # Generate some random data.
#' means <- runif(100,min=0,max=1000)
#' counts <- matrix(rnbinom(1000, size=1/0.01, mu=rep(means,10)), ncol=10)
#'
#' y <- varistran::vst(counts)
#'
#' # Information about the transformation
#' varistran::vst_advice(y)
#'
#' @export
vst <- function(x, method="anscombe.nb", lib.size=NULL, cpm=FALSE, dispersion=NULL, design=NULL) {
    x <- as.matrix(x)

    method.info <- vst_methods[[method]]
    is.null(method.info) && stop("Unknown method")

    # Empty matrix, do nothing
    if (nrow(x) == 0 || ncol(x) == 0) {
        return(x)
    }

    true.lib.size <- colSums(x)

    if (is.null(lib.size)) {
        # edgeR dies on zero library size
        good <- true.lib.size > 0
        norm.factors <- rep(1,length(good))
        norm.factors[good] <- edgeR::calcNormFactors(x[,good,drop=FALSE])
        
        lib.size <- true.lib.size * norm.factors
        lib.size.method <- "TMM normalization"
    } else if (all(lib.size == true.lib.size)) {
        lib.size.method <- "no adjustment"
    } else {
        lib.size.method <- "unknown"
    }

    # Avoid division by zero for empty sample
    lib.size <- pmax(lib.size, 1)

    mean.size <- mean(lib.size)

    x.norm <- t(t(x) * (mean.size / lib.size))

    if (method.info$needs.dispersion) {
        if (is.null(dispersion)) {
            dispersion <- .optimal_dispersion(x,method=method,lib.size=lib.size,design=design)
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
    attr(result, "true.lib.size") <- true.lib.size
    attr(result, "lib.size.method") <- lib.size.method
    attr(result, "cpm") <- cpm
    attr(result, "method") <- method
    if (method.info$needs.dispersion)
        attr(result, "dispersion") <- dispersion

    result
}


#' Advise how VST will transform data
#'
#' @param what Either the output of a call to vst() or the name of a VST method (see vst() help).
#'
#' @param dispersion As per vst().
#'
#' @param cpm As per vst().
#'
#' @param lib.size As per vst().
#'
#' @return A data frame giving an indication of how an average sample will be transformed.
#'
#' The column "twofold_step" shows the step from the previous to current row. With a log2 transformation this would be uniformly 1, but with a VST and small counts the step is less than 1. This therefore provides advice on how compacted the VST is near zero-count, as compared to a log2 transformation.
#'
#' Note that the results given are for an average sample. Where library sizes differ wildly, the VST may perform poorly.
#'
#' @export
vst_advice <- function(what="anscombe.nb", dispersion=NULL, cpm=FALSE, lib.size=NULL) {
    if (!is.character(what)) {
        (is.null(dispersion) && is.null(lib.size)) ||
            stop("Extra parameters not needed")

        method <- attr(what,"method")
        dispersion <- attr(what,"dispersion")
        cpm <- attr(what,"cpm")
        lib.size <- mean(attr(what,"lib.size"))
    } else {
        method <- what
    }

    method.info <- vst_methods[[method]]
    is.null(method.info) && stop("Unknown method")

    method.info$needs.dispersion && is.null(dispersion) && stop("dispersion needed")

    cpm && is.null(lib.size) && stop("lib.size needed")

    count <- cbind(c(0, 2**(0:12)))

    y <- vst(count, method=method,dispersion=dispersion,cpm=cpm,lib.size=lib.size)

    step <- rep(NA, nrow(count))
    step[3:nrow(count)] <- y[3:nrow(count),1] - y[2:(nrow(count)-1),1]

    data.frame(
        count=count[,1],
        transformed_count=y[,1],
        twofold_step=step)
}
