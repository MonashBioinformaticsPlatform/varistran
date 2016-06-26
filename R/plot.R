# Partition a vector evenly and calculate means of each partition
.bin <- function(vec, n) {
    result <- numeric(n)
    for(i in seq_len(n))
        result[i] <- mean(vec[ (floor((i-1)*length(vec)/n)+1):floor(i*length(vec)/n) ])
    result
}


.stack_labels_positive <- function(df, text_size) {
    width <- nchar(df$label) * text_size * 0.7
    df$left <- df$x - width*0.5
    df$right <- df$x + width*0.5
    df <- df[order(df$y),,drop=F]

    yoff <- df$y
    for(i in seq_len(nrow(df))) {
        obstacles <- seq_len(i-1)
        options <- yoff[obstacles] + text_size
        options <- options[options > df$y[i]]
        options <- sort(c(df$y[i], options))
        for(y in options) {
            good <- T
            for(j in obstacles)
                if (df$left[j] < df$right[i] &&
                    df$left[i] < df$right[j] &&
                    yoff[j] < y+text_size*0.999 &&
                    y < yoff[j]+text_size*0.999
                    )
                    good <- F

            if (good) {
                yoff[i] <- y
                break
            }
        }
    }

    df$yoff <- yoff
    df
}

# Stack up some labels
.stack_labels <- function(df, text_size) {
    top <- df$y>0

    top_rows <- df[top,,drop=F]
    top_rows$vjust <- -0.5
    top_rows <- .stack_labels_positive(top_rows, text_size)

    bottom_rows <- df[!top,,drop=F]
    bottom_rows$vjust <- 1.5
    bottom_rows$y <- -bottom_rows$y
    bottom_rows <- .stack_labels_positive(bottom_rows, text_size)
    bottom_rows$y <- -bottom_rows$y
    bottom_rows$yoff <- -bottom_rows$yoff

    rbind(
        top_rows,
        bottom_rows
    )
}


#' Stability plot.
#'
#' Produce a ggplot object containing a plot of residual standard deviation
#' against mean count.
#'
#' Genes are partitioned evenly into "bins" bins by average expression level.
#' Mean residual standard deviation is plotted against mean count.
#'
#' @param y Transformed counts matrix.
#' @param x Optional, original counts matrix.
#' @param design Matrix specifying a linear model with which to calculate
#' residuals.
#' @param bins Number points in the graph.
#' @return A ggplot object.
#'
#' This must be print()-ed to actually plot.
#' @author Paul Harrison
#'
#' @export
plot_stability <- function(y, x=NULL, design=NULL, bins=20) {
    y <- as.matrix(y)

    if (!is.null(x)) {
        x <- as.matrix(x)
        x_data <- rowMeans(x)
        x_label <- "mean count"
    } else {
        x_data <- rank(rowSums(y), ties.method="first")
        x_label <- "rank by mean expression"
    }

    if (is.null(design))
        design <- matrix(1, ncol=1, nrow=ncol(y))

    residuals <- y %*% MASS::Null(design)
    rsd <- sqrt(rowSums(residuals*residuals))

    reorder <- order(x_data)

    df <- data.frame(
        rsd=.bin(rsd[reorder], bins),
        x=.bin(x_data[reorder], bins)
    )

    result <- ggplot2::ggplot(df, ggplot2::aes_string(x="x",y="rsd")) +
        ggplot2::geom_point(size=3) + ggplot2::geom_line() +
        ggplot2::ylim(0,NA) +
        ggplot2::xlab(x_label) +
        ggplot2::ylab(
            if (ncol(design) > 1) "residual standard deviation"
            else "standard deviation"
        ) +
        ggplot2::theme_bw()

    if (!is.null(x)) {
        to <- floor(log10(max(df$x))) + 1
        result <- result + ggplot2::scale_x_log10(breaks=10^(0:to))
    }

    result
}




#' Biplot of expression data
#'
#' Produce a ggplot object containing a biplot of expression data.
#'
#' Biplot based on the Singular Value Decomposition of the matrix x. The
#' dimensions corresponding to the two largest singular values are shown.
#'
#' Genes are shown in blue and samples in red.
#'
#' The dot product of the gene and sample vectors approximates the difference
#' from the average expression level of that gene in that sample.
#'
#' Sample points (red) are scaled to have the same variance in the two
#' dimensions. Therefore the gene points (blue) may have greater variance along
#' dimension 1 if dimension 1 explains more of the variance than dimension 2.
#'
#' @param x Matrix of expression levels, with features (eg genes) as rows and
#' samples as columns. For example, you could use the output of varistran::vst
#' here.
#' @param sample_labels Sample labels.
#' @param feature_labels Feature labels.
#' @param balance Relative scaling of features and samples.
#' @param n_features Number of extreme features to label.
#' @param text_size plot_biplot attempts to stop labels from overlapping.
#' Adjust this so that text just doesn't overlap. Set to zero to allow labels
#' to completely overlap.
#' @return A ggplot object.
#'
#' This must be print()-ed to actually plot.
#' @author Paul Harrison
#' @examples
#'
#'
#' # Assuming counts is a matrix of read counts.
#' y <- varistran::vst(counts)
#' print( varistran::plot_biplot(y) )
#'
#' @export
plot_biplot <- function(x, sample_labels=NULL, feature_labels=NULL, n_features=20, balance=0.25, text_size=0.025) {
    x <- as.matrix(x)

    if (is.null(sample_labels) && !is.null(colnames(x)))
        sample_labels <- colnames(x)

    if (is.null(sample_labels))
        sample_labels <- rep("", ncol(x))

    sample_labels[is.na(sample_labels)] <- ""


    if (is.null(feature_labels) && !is.null(rownames(x)))
        feature_labels <- rownames(x)

    if (is.null(feature_labels))
        feature_labels <- rep("", nrow(x))

    feature_labels[is.na(feature_labels)] <- ""


    n_features <- min(n_features, nrow(x))

    decomp <- svd(x - rowMeans(x))

    d2 <- decomp$d ^ 2
    R2 <- d2 / sum(d2)

    balancer <- sqrt(
        sqrt(nrow(x) / ncol(x))
        / sqrt(d2[1]+d2[2])
        * balance
    )

    u <- t(t(decomp$u) * (decomp$d*balancer))
    v <- t(t(decomp$v) / balancer)

    features <- data.frame(
        x = u[,1],
        y = u[,2],
        is_feature = rep(TRUE, nrow(x)),
        label = feature_labels,
        stringsAsFactors = FALSE
    )

    samples <- data.frame(
        x = v[,1],
        y = v[,2],
        is_feature = rep(FALSE, ncol(x)),
        label = sample_labels,
        stringsAsFactors = FALSE
    )

    result <- ggplot2::ggplot(features, ggplot2::aes_string(x="x",y="y")) +
        ggplot2::coord_fixed() +
        ggplot2::geom_point(size=1.5, color="#0088ff") +
        ggplot2::geom_point(size=3,data=samples,color="red") +
        ggplot2::xlab(sprintf("Dimension 1, %.1f%% of variance", R2[1]*100)) +
        ggplot2::ylab(sprintf("Dimension 2, %.1f%% of variance", R2[2]*100)) +
        ggplot2::theme_bw()

    to_label <- samples

    score <- decomp$u[,1]^2 + decomp$u[,2]^2
    selection <- order(score,decreasing=T)[seq_len(n_features)]
    to_label <- rbind(to_label, features[selection,])

    to_label <- to_label[to_label$label != "", ]

    ylow <- min(features$y,samples$y)
    yhigh <- max(features$y,samples$y)
    xlow <- min(features$x,samples$x)
    xhigh <- max(features$x,samples$x)

    if (nrow(to_label) > 0) {
        scale <- max(
            max(features$x,samples$x)-min(features$x,samples$x),
            max(features$y,samples$y)-min(features$y,samples$y)
            )

        scaled_text_size <- text_size * scale

        to_label <- .stack_labels(to_label, scaled_text_size)

        result <- result +
            ggplot2::geom_segment(
                data=to_label,
                ggplot2::aes_string(x="x",y="y",xend="x",yend="yoff"),
                alpha=0.2)

        if (any(to_label$is_feature))
            result <- result +
                ggplot2::geom_text(
                    data=to_label[to_label$is_feature,],
                    ggplot2::aes_string(label="label",y="yoff",vjust="vjust"),
                    size=4,alpha=1/3)

        if (any(!to_label$is_feature))
            result <- result +
                ggplot2::geom_text(
                    data=to_label[!to_label$is_feature,],
                    ggplot2::aes_string(label="label",y="yoff",vjust="vjust"),
                    size=4)

        ylow <- min(ylow,min(to_label$yoff)-scaled_text_size)
        yhigh <- max(yhigh,max(to_label$yoff)+scaled_text_size)
        xlow <- min(xlow,min(to_label$left))
        xhigh <- max(xhigh,max(to_label$right))
        result <- result + ggplot2::xlim(xlow,xhigh) + ggplot2::ylim(ylow,yhigh)
    }

    result
}
