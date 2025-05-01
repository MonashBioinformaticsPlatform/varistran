
#
# Find scale to achieve log1p(x*scale) == target
#
get_scale <- function(x, target, tol=1e-9, verbose=FALSE) {
    x <- x[x!=0]
    
    if (length(x) == 0) 
        return(initial)
    
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

#' Normalized log2 counts using log2(count[i,j]*scale[j]+1)
#'
#'
#' @export
samesum_log2_norm <- function(x, scale=0.25, tol=1e-9, verbose=FALSE) {
    # Calculations are all do using natural log
    baseline <- colSums(log1p(x*scale))
    
    target <- mean(baseline)
    
    scales <- rep(0, ncol(x))
    for(i in seq_len(ncol(x))) {
        scales[i] <- get_scale(x[,i], target=target, tol=tol, verbose=verbose)
    }
    
    result <- log1p(t(t(x)*scales)) / log(2)
    
    names(scales) <- colnames(x)
    attr(result,"scale") <- scales
    
    result
}