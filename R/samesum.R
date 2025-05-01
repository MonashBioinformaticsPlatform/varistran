
#
# Find scale to achieve log1p(x*scale) == target
#
get_scale <- function(x, target, tol=1e-9, verbose=FALSE) {
    x <- x[x!=0]
    
    if (length(x) == 0) 
        return(0)
    
    scale <- 0
    iter <- 0
    repeat {
        value <- sum(log1p(x*scale))
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
        deriv <- sum(x/(x*scale+1))
        scale <- scale - (value-target)/deriv
    }
    
    if (verbose) {
        message("Iterations: ", iter)
    }
    
    scale
}

#' Normalized log2 counts using "scalesum" method
#'
#' @param scale Aim to give each sample close to this value of scale.
#'
#' @param target_sum Aim to give each sample exactly this sum. "scale" is ignored if this parameter is used. This could be used to provide a transformation that is consistent across datasets.
#'
#' @export
samesum_log2_norm <- function(x, scale=0.25, target_sum=NA, tol=1e-9, verbose=FALSE) {
    # Calculations are all done using natural log
    
    if (is.na(target_sum)) {
        baseline <- colSums(log1p(x*scale))
        target_sum <- mean(baseline)
    } else {
        target_sum <- target_sum * log(2)
    }
    
    scales <- rep(0, ncol(x))
    for(i in seq_len(ncol(x))) {
        scales[i] <- get_scale(x[,i], target=target_sum, tol=tol, verbose=verbose)
    }
    
    # Comment on anything concerning
    
    n_bad <- sum(scales >= 1 & scales != Inf)
    if (n_bad) {
        warning(n_bad, " samples have been given a scale greater than one. If these subsequently appear to be outliers, this may be an artifact of the log1p transformation. Consider using a smaller target scale or excluding these samples from analysis.")
    }
    
    n_inf <- sum(scales == Inf)
    if (n_inf) {
        warning(n_inf, " samples where the scale has overflowed and been given the value Inf. Please remove these samples.")
    }
    
    sums <- colSums(x)
    n_zero <- sum(sums == 0)
    if (n_zero > 0) {
        warning(n_zero, " samples are entirely zero.")
    }
    
    
    # Assemble result
    
    result <- log1p(t(t(x)*scales)) / log(2)
    
    samples <- data.frame(row.names=colnames(x) %||% seq_len(ncol(x)))
    samples$name <- colnames(x)
    samples$scale <- scales
    
    scaled_sums <- sums * samples$scale
    norm_factor <- scaled_sums/exp(mean(log(scaled_sums))) 
    samples$true_lib_size <- sums
    samples$adjusted_lib_size <- sums * norm_factor 
    samples$norm_factor <- norm_factor
    
    attr(result,"samples") <- samples
    attr(result,"zero") <- 0
    
    result
}

#' Normalized log2 CPM using "scalesum" method
#'
#' Applies samesum_log2_norm, then adds an adjustment to produce values that can be treated as log2 Counts Per Million (CPM).
#'
#' The result will have the property that \code{mean(colSums(2**lcpm-2**attr(lcpm,"zero")))} equals one million.
#'
#' @export
samesum_log2_cpm <- function(x, scale=0.25, target_sum=NA, tol=1e-9, verbose=FALSE) {
    result <- samesum_log2_norm(x, scale=scale, target_sum=target_sum, tol=tol, verbose=verbose)
    
    samples <- attr(result,"samples")
    scaled_sums <- samples$true_lib_size * samples$scale
    adj <- log2(1e6) - log2(mean(scaled_sums)) 
    
    result <- result + adj
    attr(result, "zero") <- adj
    
    result
}