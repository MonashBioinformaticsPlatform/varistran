

plot_biplot <- function(x, sample_labels=NULL, balance=0.25) {
    x <- as.matrix(x)
    
    if (is.null(sample_labels) && !is.null(colnames(x)))
        sample_labels <- colnames(x)
    
    decomp <- svd(x - rowMeans(x))

    balancer <- ( nrow(x) / ncol(x) )**0.25 * sqrt(balance)

    u <- t(t(decomp$u) * (sqrt(decomp$d)*balancer)) 
    v <- t(t(decomp$v) * (sqrt(decomp$d)/balancer)) 
    
    d2 <- decomp$d ^ 2
    R2 <- d2 / sum(d2)

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
            geom_text(data=samples,aes(label=labels), vjust=ifelse(samples$y<0,1.25,-0.25))
        
    result <- result +
        geom_point(size=4,data=samples,color="red") +
        xlab(sprintf("Dimension 1, %.1f%% of variance", R2[1]*100)) +
        ylab(sprintf("Dimension 2, %.1f%% of variance", R2[2]*100))

    result
}

