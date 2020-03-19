
#' Plot a heatmap.
#'
#' Produces a heatmap as a grid grob.
#'
#' This heatmap differs from other heatmaps in R in the method of clustering used:
#'
#' 1. The distances used are cosine distances (i.e. the magnitude of log fold changes is not important, only the pattern).
#'
#' 2. \code{hclust()} is used to produce a clustering, as normal.
#'
#' 3. Branches in the hierarchical clustering are flipped to minimize sharp changes between neighbours, using the seriation package's OLO (Optimal Leaf Ordering) method.
#'
#' @param y A matrix of expression levels, such as a transformed counts matrix as produced by \code{varistran::vst}.
#' @param cluster_samples Should samples (columns) be clustered?
#' @param cluster_features Should features (rows) be clustered?
#' @param sample_labels Names for each sample. If not given and y has column names, these will be used instead.
#' @param feature_labels Names for each feature. If not given and y has row names, these will be used instead.
#' @param baseline Baseline level for each row, to be subtracted when drawing the heatmap colors. If omitted, the row mean will be used. Specify \code{baseline=0} to not subtract anything and not show a baseline bar graph.
#' @param baseline_label Text description of what the baseline is.
#' @param scale_label Text description of what the heatmap colors represent (after baseline is subtracted).
#' @param n Show only this many rows. Rows are selected in order of greatest span of expression level.
#'
#' @return A grid grob. print()-ing this value will cause it to be displayed.
#'
#' Additionally $info$row_order will contain row ordering and $info$col_order will contain column ordering.
#' @author Paul Harrison.
#'
#' @examples
#'
#' # Generate some random data.
#' counts <- matrix(rnbinom(1000, size=1/0.01, mu=100), ncol=10)
#'
#' y <- varistran::vst(counts, cpm=TRUE)
#' print( varistran::plot_heatmap(y, n=20) )
#'
#'
#' @export
plot_heatmap <- function(
        y,
        cluster_samples=FALSE,
        cluster_features=TRUE,
        sample_labels=NULL,
        feature_labels=NULL,
        baseline=NULL,
        baseline_label="row\nmean",
        scale_label="difference from\nrow mean",
        n=Inf) {
    y <- as.matrix(y)

    if (is.null(sample_labels) && !is.null(colnames(y)))
        sample_labels <- colnames(y)

    if (is.null(sample_labels))
        sample_labels <- rep("", ncol(y))

    sample_labels[is.na(sample_labels)] <- ""


    if (is.null(feature_labels) && !is.null(rownames(y)))
        feature_labels <- rownames(y)

    if (is.null(feature_labels))
        feature_labels <- rep("", nrow(y))

    feature_labels[is.na(feature_labels)] <- ""

    if (!is.null(baseline)) {
        if (length(baseline) == 1)
            baseline <- rep(baseline, nrow(y))
        stopifnot(length(baseline) == nrow(y))
        means <- baseline
    } else {
        means <- rowMeans(y, na.rm=TRUE)
    }

    # Show only a subset of rows, if desired
    if (n < nrow(y)) {    
        y_span <- apply(y,1,max,-Inf,na.rm=TRUE) - apply(y,1,min,Inf,na.rm=TRUE)
        y_span[ !is.finite(y_span) ] <- -Inf
        selection <- rep(FALSE,nrow(y))
        selection[ order(-y_span)[ seq_len(n) ] ] <- TRUE

        y <- y[selection,,drop=FALSE]
        feature_labels <- feature_labels[selection]
        means <- means[selection]
    }


    y_centered <- y - means
    
    y_scaled <- y_centered / sqrt(rowMeans(y_centered*y_centered, na.rm=TRUE))
    y_scaled[ is.na(y_scaled) ] <- 0.0
    
    row_order <- make_ordering(y_scaled, enable=cluster_features)

    y_centered_clean <- y_centered
    y_centered_clean[ is.na(y_centered_clean) ] <- 0.0
    col_order <- make_ordering(t(y_centered_clean), enable=cluster_samples)

    pad <- 0.25

    row_ordering_grob <- ordering_grob(row_order, transpose=TRUE, mirror=TRUE)

    col_ordering_grob <- ordering_grob(col_order)

    heatmap <- heatmap_grob(
        y_centered[row_order$order,col_order$order,drop=F],
        signed=TRUE,
        legend_title=paste0(scale_label),
        vp_name="heatmap")

    mean_range <- range(means, na.rm=TRUE)
    
    need_means <- mean_range[1] != 0 || mean_range[2] != 0
    
    if (mean_range[2] == mean_range[1]) 
        mean_range[2] <- mean_range[2]+1
    
    if (need_means) {
        mean_graph <- rectGrob(
            x=rep(mean_range[1],nrow(y)),
            y=seq_len(nrow(y))-1,
            width=means[row_order$order]-mean_range[1],
            height=rep(1,nrow(y)),
            just=c(0,0),
            default.units="native",
            vp=viewport(xscale=mean_range,yscale=c(0,nrow(y)))
        )
        mean_axis <- xaxisGrob(
            at=axisTicks(mean_range,log=FALSE,nint=3),
            label=TRUE,
            vp=viewport(width=1,height=0,y=1,xscale=mean_range),
            gp=gpar(cex=0.75)
        )
        mean_label <- textGrob(baseline_label)
        mean_width <- unit(3,"lines")
        mean_pad <- pad
    } else {
        mean_graph <- textGrob("")
        mean_axis <- textGrob("")
        mean_label <- textGrob("")
        mean_width <- unit(0,"lines")
        mean_pad <- 0
    }
    
    feature_label_grob <- shrinktext_grob(
        feature_labels[row_order$order],
        x=rep(0,nrow(y)),
        y=seq_len(nrow(y))-0.5,
        just=c(0,0.5),
        vp=viewport(xscale=c(0,1),yscale=c(0,nrow(y)))
    )

    sample_label_grob <- vertical_shrinktext_grob(
        sample_labels[col_order$order],
        x=seq_len(ncol(y))-0.5,
        y=rep(1,ncol(y)),
        just=c(1,0.5),
        vp=viewport(xscale=c(0,ncol(y)),yscale=c(0,1))
    )

    frame <- frameGrob(layout=grid.layout(nrow=3,ncol=4))

    frame <- packGrob(frame, varistran_grob(col_ordering_grob,height="inherit",pad=pad), row=1,col=2)
    frame <- packGrob(frame, varistran_grob(mean_label,height="inherit",pad=mean_pad), row=1,col=3)

    frame <- packGrob(frame, varistran_grob(row_ordering_grob,width="inherit",pad=pad), row=2,col=1)
    frame <- packGrob(frame, varistran_grob(heatmap$heatmap,pad=pad), row=2, col=2)
    frame <- packGrob(frame, varistran_grob(mean_graph,width=mean_width,pad=mean_pad), row=2,col=3)
    frame <- packGrob(frame, varistran_grob(feature_label_grob,width="inherit",pad=pad), row=2,col=4)

    frame <- packGrob(frame, varistran_grob(sample_label_grob,height="inherit",pad=pad), row=3,col=2)
    frame <- packGrob(frame, varistran_grob(mean_axis,height=unit(3,"lines"),pad=mean_pad), row=3,col=3)

    outer <- frameGrob()
    outer <- packGrob(outer, varistran_grob(frame), row=1,col=1)
    outer <- packGrob(outer, varistran_grob(heatmap$legend,height="inherit",pad=pad), row=2, col=1)

    result <- varistran_grob(outer, pad=pad)
    result$info <- list(
        row_order=row_order,
        col_order=col_order
    )
    result
}
