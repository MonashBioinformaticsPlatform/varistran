
# asinh(ab) / log(2) + const, behaving as log2(a) in the limit as a -> large
.log2ish.asinh <- function(a,b) {
    ab <- a * b
    log2(ab+sqrt(1+ab*ab)) - log2(b) - 1
}


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
    
    residuals <- mat %*% MASS::Null(design)
    
    rss <- rowSums(residuals*residuals)
    
    sqrt(mean(rss)) / mean(sqrt(rss))
}


.optimal.dispersion <- function(x, method="anscombe.nb", lib.size=NULL, design=NULL) {
    x <- x[ rowMeans(x) >= 5, ]

    optimize(
         function(d) {             
             .dispersion.score(
                 vst(x,method=method,lib.size=NULL,dispersion=d), 
                 design=design)
         },
         lower = 1e-4,
         upper = 1.0
    )$minimum
}


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
    
    attr(result, "lib.size") <- lib.size
    attr(result, "cpm") <- cpm
    attr(result, "method") <- method
    if (method.info$needs.dispersion)
        attr(result, "dispersion") <- dispersion

    result
}






