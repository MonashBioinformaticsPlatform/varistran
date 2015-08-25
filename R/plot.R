

# Partition a vector evenly and calculate means of each partition
.bin <- function(vec, n) {
    result <- numeric(n)
    for(i in seq_len(n))
        result[i] <- mean(vec[ (floor((i-1)*length(vec)/n)+1):floor(i*length(vec)/n) ])
    result
}


# Stack some labels
.stack_labels <- function(df, text_size) {
    df$vjust <- ifelse(df$y<0,1.5,-0.5)
    width <- nchar(df$label) * text_size * 0.75
    df$left <- df$x - width*0.5
    df$right <- df$x + width*0.5
    df$yoff <- df$y
    
    stack_order <- order(abs(df$y))
    for(i in seq_along(stack_order)) {
       soi <- stack_order[i]
       for(j in seq_len(i-1)) {
           soj <- stack_order[j]
           if (df$vjust[soi] == df$vjust[soj] && 
               df$left[soj] < df$right[soi] && 
               df$left[soi] < df$right[soj]) {
               if (df$y[soi]<0) {
                   df$yoff[soi] <- min(df$yoff[soi], df$yoff[soj]-text_size) 
               } else {
                   df$yoff[soi] <- max(df$yoff[soi], df$yoff[soj]+text_size) 
               }
           }
       }
    }
    
    df
} 



plot_stability <- function(y, x=NULL, design=NULL, bins=20) {
    y <- as.matrix(y)
    
    if (!is.null(x)) {
        x <- as.matrix(x)
        x_data <- rowMeans(x)
        x_label <- "mean count"
    } else {
        x_data <- rank(rowSums(y), ties="first")
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
    
    result <- ggplot2::ggplot(df, ggplot2::aes(x=x,y=rsd)) +
        ggplot2::geom_point(size=3) + ggplot2::geom_line() +
        ggplot2::ylim(0,NA) +
        ggplot2::xlab(x_label) +
        ggplot2::ylab("residual standard deviation") +
        ggplot2::theme_bw()
   
    if (!is.null(x)) {
        to <- floor(log10(max(df$x))) + 1
        result <- result + ggplot2::scale_x_log10(breaks=10^(0:to))
    }
    
    result
}


plot_biplot <- function(x, sample_labels=NULL, feature_labels=NULL, n_features=20, balance=0.25, text_size=0.025) {
    x <- as.matrix(x)
    
    if (is.null(sample_labels) && !is.null(colnames(x)))
        sample_labels <- colnames(x)

    if (is.null(feature_labels) && !is.null(rownames(x)))
        feature_labels <- rownames(x)

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
        is_feature = rep(TRUE, nrow(x))
    )
    if (!is.null(feature_labels))
        features$labels <- feature_labels
    
    samples <- data.frame(
        x = v[,1],
        y = v[,2],
        is_feature = rep(FALSE, ncol(x))
    )
    if (!is.null(sample_labels)) {
        samples$labels <- sample_labels
        to_label <- samples
    } else {
        to_label <- NULL
    }
    
    result <- ggplot2::ggplot(features, ggplot2::aes(x=x,y=y)) + 
        ggplot2::coord_fixed() + 
        ggplot2::geom_point(size=1.5, color="#0088ff") +
        ggplot2::geom_point(size=3,data=samples,color="red") +
        ggplot2::xlab(sprintf("Dimension 1, %.1f%% of variance", R2[1]*100)) +
        ggplot2::ylab(sprintf("Dimension 2, %.1f%% of variance", R2[2]*100)) +
        ggplot2::theme_bw()

    if (!is.null(feature_labels) && n_features > 0) {
        score <- decomp$u[,1]^2 + decomp$u[,2]^2
        selection <- order(score,decreasing=T)[seq_len(n_features)]
        
        if (is.null(to_label))
            to_label <- features[selection,]
        else
            to_label <- rbind(to_label, features[selection,])
    }
    
    ylow <- min(features$y,samples$y)
    yhigh <- max(features$y,samples$y)
    xlow <- min(features$x,samples$x)
    xhigh <- max(features$x,samples$x)
    
    if (!is.null(to_label)) {
        scale <- max(
            max(features$x,samples$x)-min(features$x,samples$x),
            max(features$y,samples$y)-min(features$y,samples$y)
            )
        
        scaled_text_size <- text_size * scale
        
        to_label <- .stack_labels(to_label, scaled_text_size)
        
        result <- result +
            ggplot2::geom_segment(data=to_label, ggplot2::aes(x=x,y=y,xend=x,yend=yoff), alpha=0.2)

        #result <- result +
        #    ggplot2::geom_segment(data=to_label, ggplot2::aes(x=left,y=yoff,xend=right,yend=yoff), color="#888888")
            
        if (any(to_label$is_feature))
            result <- result + 
                ggplot2::geom_text(data=to_label[to_label$is_feature,],ggplot2::aes(label=labels,y=yoff,vjust=vjust), alpha=0.5)
        
        if (any(!to_label$is_feature))
            result <- result + 
                ggplot2::geom_text(data=to_label[!to_label$is_feature,],ggplot2::aes(label=labels,y=yoff,vjust=vjust))
        
        ylow <- min(ylow,min(to_label$yoff)-scaled_text_size)
        yhigh <- max(yhigh,max(to_label$yoff)+scaled_text_size)
        xlow <- min(xlow,min(to_label$left))
        xhigh <- max(xhigh,max(to_label$right)) 
        result <- result + ggplot2::xlim(xlow,xhigh) + ggplot2::ylim(ylow,yhigh)
    }
    
    result
}

