
#' @export
plot_heatmap <- function(
        y,
        cluster_samples=FALSE,
        cluster_features=TRUE,
        sample_labels=NULL,
        feature_labels=NULL,
        units="units") {
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


    means <- rowMeans(y)
    y_centered <- y - means

    y_scaled <- y_centered / sqrt(rowMeans(y_centered*y_centered))
    row_order <- make_ordering(y_scaled, enable=cluster_features)

    col_order <- make_ordering(t(y_centered), enable=cluster_samples)

    pad <- 0.25

    row_ordering_grob <- ordering_grob(row_order, transpose=TRUE, mirror=TRUE)

    col_ordering_grob <- ordering_grob(col_order)

    heatmap <- heatmap_grob(y_centered[row_order$order,col_order$order,drop=F], signed=TRUE, legend_title=paste0(units,",\ndifference from row mean"))

    mean_range <- range(means)
    if (mean_range[2] == mean_range[1]) mean_range[2] <- mean_range[2]+1
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
    mean_label <- textGrob("row\nmean")

    feature_label_grob <- textGrob(
        feature_labels[row_order$order],
        x=rep(0,nrow(y)),
        y=seq_len(nrow(y))-0.5,
        just=c(0,0.5),
        default.units="native",
        vp=viewport(xscale=c(0,1),yscale=c(0,nrow(y)))
    )

    sample_label_grob <- textGrob(
        sample_labels[col_order$order],
        x=seq_len(ncol(y))-0.5,
        y=rep(1,ncol(y)),
        just=c(1,0.5),
        rot=90,
        default.units="native",
        vp=viewport(xscale=c(0,ncol(y)),yscale=c(0,1))
    )

    frame <- frameGrob(layout=grid.layout(nrow=3,ncol=4))

    frame <- packGrob(frame, varistran_grob(col_ordering_grob,height="inherit",pad=pad), row=1,col=2)
    frame <- packGrob(frame, varistran_grob(mean_label,height="inherit",pad=pad), row=1,col=3)

    frame <- packGrob(frame, varistran_grob(row_ordering_grob,width="inherit",pad=pad), row=2,col=1)
    frame <- packGrob(frame, varistran_grob(heatmap$heatmap,pad=pad), row=2, col=2)
    frame <- packGrob(frame, varistran_grob(mean_graph,width=unit(3,"lines"),pad=pad), row=2,col=3)
    frame <- packGrob(frame, varistran_grob(feature_label_grob,width="inherit",pad=pad), row=2,col=4)

    frame <- packGrob(frame, varistran_grob(sample_label_grob,height="inherit",pad=pad), row=3,col=2)
    frame <- packGrob(frame, varistran_grob(mean_axis,height=unit(3,"lines"),pad=pad), row=3,col=3)

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
