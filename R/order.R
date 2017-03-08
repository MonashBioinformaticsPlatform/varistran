
#
# Represent an ordering, primarily for heatmaps
#

dendrogram_paths <- function(dend) {
   if (is.leaf(dend)) {
       ''
   } else {
       result <- character(0)
       for(i in 1:length(dend)) {
           result <- c(result, paste(sprintf('%d',length(dend)-i), dendrogram_paths(dend[[i]])))
       }
       result
   }
}

#' Make an object representing an ordering, usually from a hierarchical clustering.
#' 
#' Performs hierarchical clustering, and then further refines the order using the seriation library.
#'
#' @param mat Matrix suitably transformed so distances between rows are meaningful.
#'
#' @param enable If false, don't perform clustering.
#'
#' @param fast If true, just order rows with a simple one dimensional PCA.
#'
#' @export
make_ordering <- function(mat, enable=TRUE, fast=FALSE) {
    # Note: paths are given ordered by order

    # Fill in NAs
    fillin <- !is.finite(mat)
    if (any(fillin)) {
        filler <- matrix(rowMeans(mat, na.rm=TRUE), nrow=nrow(mat),ncol=ncol(mat))
        mat[fillin] <- filler[fillin]
        mat[!is.finite(mat)] <- mean(mat)
        mat[!is.finite(mat)] <- 0.0
    }
    
    if (nrow(mat) < 3 || !enable) {
        list(
            dendrogram = NULL,
            order = seq_len(nrow(mat)),
            paths = rep('',nrow(mat)))

    } else if (fast) {
        # Too many for seriation to deal with efficiently
        perm <- seriation::seriate(mat, method='PCA')[[1]]

        list(dendrogram = NULL,
             order = seriation::get_order(perm),
             paths = rep('',nrow(mat)))

    } else {
        dist_mat <- dist(mat)
        control <- list(hclust = hclust(dist_mat))
        dend_mat <- as.dendrogram(
                seriation::seriate(dist_mat,
                method = 'OLO',
                control = control)[[1]])

        list(dendrogram = dend_mat,
             order = order.dendrogram(dend_mat),
             paths = dendrogram_paths(dend_mat))
    }
}

