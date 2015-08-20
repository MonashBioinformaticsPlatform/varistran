

# Partition a vector evenly and calculate means of each partition
.bin <- function(vec, n) {
    result <- numeric(n)
    for(i in seq_len(n))
        result[i] <- mean(vec[ (floor((i-1)*length(vec)/n)+1):floor(i*length(vec)/n) ])
    result
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
    print(df)
    
    result <- ggplot(df, aes(x=x,y=rsd)) +
        geom_point(size=3) + geom_line() +
        xlab(x_label) +
        ylab("residual standard deviation")
   
    if (!is.null(x)) {
        to <- floor(log10(max(df$x))) + 1
        result <- result + scale_x_log10(breaks=10^(0:to))
    }
    
    result
}


plot_biplot <- function(x, sample_labels=NULL, balance=0.25) {
    x <- as.matrix(x)
    
    if (is.null(sample_labels) && !is.null(colnames(x)))
        sample_labels <- colnames(x)
    
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
        y = u[,2]
    )
    
    samples <- data.frame(
        x = v[,1],
        y = v[,2]
    )
    if (!is.null(sample_labels))
        samples$labels <- sample_labels
    
    
    result <- ggplot(features, aes(x=x,y=y)) + 
        coord_fixed() + 
        geom_point(size=1, color="#0088ff")
        
    if (!is.null(sample_labels))
        result <- result + 
            geom_text(data=samples,aes(label=labels), vjust=ifelse(samples$y<0,1.5,-0.5))
        
    result <- result +
        geom_point(size=4,data=samples,color="red") +
        xlab(sprintf("Dimension 1, %.1f%% of variance", R2[1]*100)) +
        ylab(sprintf("Dimension 2, %.1f%% of variance", R2[2]*100))

    result
}

