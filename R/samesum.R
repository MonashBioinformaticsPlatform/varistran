
#
# Find scale to achieve log1p(x/scale) == target
#
get_scale <- function(x, target, tol=1e-9, verbose=FALSE) {
    x <- x[x!=0]
    
    if (length(x) == 0) 
        return(Inf)
    
    iscale <- 0
    iter <- 0
    repeat {
        value <- sum(log1p(x*iscale))
        if (value >= target*(1-tol)) {
            break
        }
        
        iter <- iter+1
        if (iter >= 100) {
            # This should never happen.
            # Newton's method should be very very fast.
            warning("Scale optimization failed to converge after 100 iterations.")
            break
        }
        
        # Newton's method step
        deriv <- sum(x/(x*iscale+1))
        iscale <- iscale - (value-target)/deriv
    }
    
    if (verbose) {
        message("Iterations: ", iter)
    }
    
    1/iscale
}

#' Normalized log2 counts using samesum method
#'
#' This is a method of computing log counts while dealing sensibly with zeros and differing library sizes. It should cope well with very sparse data. The input is transformed using \code{log2(x/scale+1)} with a different scale for each sample. The scale for each sample is chosen so all of the samples add up to the same value.
#'
#' You can either specify the argument \code{scale}, and it will try to give individual samples scale values close to this value, or directly specify the target sum with \code{sum}.
#'
#' \code{samesum_log2_cpm} adds a constant to the result to produce values that can be treated as log2 Counts Per Million (CPM). The result will have the property that \code{mean(colSums(2^lcpm-2^attr(lcpm,"zero")))} equals one million.
#'
#' Unlike using library size normalization, the samesum method is robust to variations in high abundance features, because it works on the log transformed values. It is also robust to noise from low abundance features because of the \code{+1} in the transformation.
#'
#' This method has some similarity to Centered Log Ratio (CLR) transformation. Where the CLR ensures each sample adds to zero, the samesum transformation ensures each sample adds to the same non-zero value. Like the CLR, distances between samesum transformed data are a good way to measure sample similarity. Transformed values could be used with PCA or with non-linear dimension reduction methods.
#'
#' Samesum transformed data can also be used in differential abundance analysis using the limma-trend or limma-vooma methods.
#'
#' The default value of 2 for scale should be reasonable for count data, but is not variance stabilizing. limma-trend or limma-vooma will take this into account. For visualization and exploratory methods you could try a larger scale value, which will damp down variation in low abundance rows and provide better variance stabilization. The scale parameter acts similarly to a pseudocount in log transformation methods with a pseudocount parameter (such as edgeR's \code{cpm} function).
#'
#' When different samples have different library sizes, some of them may be given a scale much smaller than the target scale. When this happens, these samples may appear as outliers in data visualization, purely as an artifact of the transformation. A warning is produced if any samples are given a scale less than one. Consider excluding these samples from your analysis. If all the samples need to be used, you can increase the target scale to avoid this problem. This will allow you to compare all samples, but at the cost of losing some ability to resolve fine differences between samples.
#'
#' @param x A matrix or a sparse matrix. Usually these will be counts, but any matrix of non-negative values could potentially be used. We regard columns as samples.
#'
#' @param scale Target scale. Aim to give each sample close to this value of scale. The default value of 2 should be reasonable for count data, but see discussion below.
#'
#' @param sum Target sum. Aim to give each sample exactly this sum. \code{scale} is ignored if this parameter is used. This could be used to provide a transformation that is consistent across datasets.
#'
#' @param tol Tolerance when optimizing scale values. Scale is increased using Newton's method until the sum is at least \code{target_sum*(1-tol)}.
#'
#' @param verbose If TRUE, output some debugging messages.
#'
#' @return 
#' A matrix or a sparse matrix, depending on the input.
#'
#' The return value will have some attributes accessible using \code{attr()}. 
#'
#' \itemize{
#' \item "samples" is a data frame with per-sample information. The most important value for the transformation is the "scale" column. There is also "true_library_size", which is the sum of the input column, and "adjusted_library_size" if you want to use this method just as a library size adjustment method. 
#' \item "zero" is the value that zero is transformed to. This will be zero for \code{samesum_log2_norm}.
#' \item "sum" is the value each column of the output sums to (but not including the CPM adjustment if using \code{samesum_log2_cpm}).
#' }
#'
#' @export
samesum_log2_norm <- function(x, scale=2, sum=NA, tol=1e-9, verbose=FALSE) {
    # Calculations are all done using natural log
    
    if (is.na(sum)) {
        baseline <- colSums(log1p(x/scale))
        sum <- mean(baseline)
    } else {
        sum <- sum * log(2)
    }
    
    scales <- rep(0, ncol(x))
    for(i in seq_len(ncol(x))) {
        scales[i] <- get_scale(x[,i], target=sum, tol=tol, verbose=verbose)
    }
    
    # Comment on anything concerning
    
    n_bad <- sum(scales < 1 & scales != 0)
    if (n_bad) {
        warning(n_bad, " samples have been given a scale less than one. If these subsequently appear to be outliers, this may be an artifact of the log1p transformation. Consider using a smaller target scale or excluding these samples from analysis.")
    }
    
    n_zero <- sum(scales == 0)
    if (n_zero) {
        warning(n_zero, " samples where the scale has underflowed and been given a value of zero. Please remove these samples.")
    }
    
    sums <- colSums(x)
    n_zero_lib <- sum(sums == 0)
    if (n_zero_lib > 0) {
        warning(n_zero_lib, " samples are entirely zero.")
    }
    
    
    # Assemble result
    
    result <- log1p(t(t(x)/scales)) / log(2)
    
    samples <- data.frame(row.names=seq_len(ncol(x)))
    samples$name <- colnames(x)
    samples$scale <- scales
    
    scaled_sums <- sums / samples$scale
    good <- is.finite(scaled_sums) & scaled_sums != 0
    norm_factor <- scaled_sums/exp(mean(log(scaled_sums[good]))) 
    samples$true_lib_size <- sums
    samples$adjusted_lib_size <- sums * norm_factor 
    
    # Thought about including this to facilite comparison to calcNormFactors, but I think it's too confusing.
    #samples$norm_factor <- norm_factor
    
    attr(result,"samples") <- samples
    attr(result,"zero") <- 0
    attr(result,"sum") <- sum / log(2)
    
    result
}

#' @rdname samesum_log2_norm
#' @export
samesum_log2_cpm <- function(x, scale=2, sum=NA, tol=1e-9, verbose=FALSE) {
    result <- samesum_log2_norm(x, scale=scale, sum=sum, tol=tol, verbose=verbose)
    
    samples <- attr(result,"samples")
    scaled_sums <- samples$true_lib_size / samples$scale
    good <- is.finite(scaled_sums) & scaled_sums != 0
    adj <- log2(1e6) - log2(mean(scaled_sums[good]))
    
    result <- result + adj
    attr(result, "zero") <- adj
    
    result
}