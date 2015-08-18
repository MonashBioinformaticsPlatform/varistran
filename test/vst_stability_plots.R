
source("test/common.R")

library("MASS")
library("Biobase")
library("edgeR")
library("DESeq2")


bin <- function(vec, n) {
    result <- numeric(n)
    for(i in seq_len(n))
        result[i] <- mean(vec[ (floor((i-1)*length(vec)/n)+1):floor(i*length(vec)/n) ])
    result
}


stability_plot <- function(prefix, x, design, dispersion, y_limit=NA, blind=F) {
    cat(prefix, nrow(counts), "genes\n")
    tnull <- t(Null(design))
    means <- rowMeans(x)
    totals <- rowSums(x)
    ord <- order(totals)    
    
    design_for_estimation <- (if (blind) NULL else design)

    ys = list(
       "log(x+c)" = 
           vst(x, method="anscombe.nb.simple", 
               dispersion=dispersion, cpm=T, design=design_for_estimation)
       ,
       "Anscombe" = 
           vst(x, method="anscombe.nb", 
               dispersion=dispersion, cpm=T, design=design_for_estimation)
       ,
       "Naive" = 
           vst(x, method="naive.nb", 
               dispersion=dispersion, cpm=T, design=design_for_estimation)
       #,
       #"cpm" = cpm(x, prior.count=0.5/dispersion, log=T)
    )
    
    dfs <- list()
    
    for(name in names(ys)) {
        #print(name)
        y <- ys[[name]]
        #print(.dispersion.score(y, design))
        variances <- apply(y,1,function(x) {
           r <- tnull %*% x
           mean(r*r)
        })
    
        
        n <- 20
        dfs[[name]] <- data.frame(
            name = rep(name,nrow(binned)),
            means = bin(means[ord], n),
            variances = bin(variances[ord], n),
            sds = bin(sqrt(variances[ord]), n)
            )
    }
    
    bigdf <- do.call(rbind, dfs)

    bigdf$name <- factor(bigdf$name, names(ys))
    
    from <- 0 #floor(log10(min(bigdf$means)))
    to <- floor(log10(max(bigdf$means))) + 1
    
    save_ggplot(
        prefix,
        ggplot(bigdf, aes(x=means,y=sds,color=name)) 
        + geom_line() 
        + geom_point(size=3) 
        + scale_x_log10(breaks=10^(from:to))
        + ylim(0.0, y_limit)
        + labs(x="mean count", y="residual standard deviation", color="Transformation")
        + my_theme
        + theme(legend.justification=c(1,0), legend.position=c(1,0))
        )
      
    #plot(rank(totals,ties.method="random"),sqrt(variances),cex=0.1)

}


real_stability_plot <- function(prefix, counts, design, y_limit=NA) {
    dgelist <- DGEList(counts)
    dgelist <- calcNormFactors(dgelist)
    dgelist <- estimateGLMCommonDisp(dgelist, design)
    
    dispersion <- dgelist$common.dispersion
    
    cat("EdgeR estimates the dispersion as",dispersion,"\n")
    
    # ylimit set so plots are comparable
    
    stability_plot(prefix, counts, design, dispersion, y_limit=y_limit)
    stability_plot(paste0(prefix,"_est"), counts, design, NULL, y_limit=y_limit)
    stability_plot(paste0(prefix,"_est_blind"), counts, design, NULL, y_limit=y_limit, blind=T)
} 

# NBPSeq arab dataset

if (T) {
    library("NBPSeq")
    data("arab")
    counts <- arab
    counts <- counts[rowSums(counts) >= 1, ]
    group <- factor(c("a","a","a","b","b","b"))
    batch <- factor(c("b1","b2","b3","b1","b2","b3"))
    design <- model.matrix(~ group + batch)
    real_stability_plot("stability_arab", counts, design, y_limit=0.65)
}


# Bottomly data set, as provided by ReCount
# Compares expression of two common inbred mouse strains

if (T) {
    bottomly.eset <- load_bottomly()
    
    counts <- exprs(bottomly.eset)
    counts <- counts[rowSums(counts) >= 1, ]
    strain <- phenoData(bottomly.eset)$strain
    experiment.number <- factor( phenoData(bottomly.eset)$experiment.number )
    design <- model.matrix(~ strain + experiment.number)

    real_stability_plot("stability_bottomly", counts, design, y_limit=0.3)
}

# Simulated data

if (T) {
    dispersion <- 0.1
    means <- 10 ^ seq(from=0,to=4,length.out=100000)
    nsamples <- 4
    #groups <- factor(c("1","1","2","2"))
    #nsamples <- length(groups)
    
    # Always produce same output
    set.seed(2015)
    
    # Matrix update in place doesn't work, so use vector.
    # library(pryr)
    counts <- numeric(length(means) * length(groups))
    for(i in seq_along(means)) {
        counts[((i-1)*nsamples+1):(i*nsamples)] <- 
            rnbinom(nsamples, size=1/dispersion, mu=means[i])
    
        # if (i%%1000 == 0) print(c(address(counts), refs(counts)))
    }
    
    counts <- matrix(counts, byrow=T,ncol=4)
    colnames(counts) <- c("a","b","c","d")
    
    #design.formula <- ~ groups
    #design <- model.matrix(design.formula)
    design <- matrix(1, ncol=1,nrow=nsamples)
    
    stability_plot("stability_simulated", counts, design, dispersion, y_limit=0.55)
    stability_plot("stability_simulated_est", counts, design, NULL, y_limit=0.55)
    stability_plot("stability_simulated_est_blind", counts, design, NULL, y_limit=0.55, blind=T)
}


#
#tnull <- t(Null(design))
#y <- DGEList(counts)
#y <- calcNormFactors(y)
#y <- estimateGLMCommonDisp(y, design)
#y <- estimateGLMTrendedDisp(y, design)
#
#dispersion <- y$common.dispersion
##dispersion <- min(y$trended.dispersion)
#print(dispersion)
#
#
##dispersion <- 1.0 #dispersion * 0.4
#print(dispersion)
#
##dispersion <- optimal.dispersion(counts, design=design)
##print(dispersion)
#
##dispersion <- NULL
#
##dds <- DESeqDataSetFromMatrix(counts, colData=data.frame(name=colnames(counts),group=group), #design = design.formula)
##dds <- DESeq(dds)
##print(dds@dispersionFunction)
#
#
##v <- vst(counts,method="anscombe.nb",dispersion=dispersion,cpm=T)
##plotMDS(v, top=nrow(v))
##
##low <- apply(v,1,min)
##high <- apply(v,1,max)
##span <- high-low
##good <- order(span, decreasing=T)[1:500]
##nesoni.heatmap(t(scale(t(v[good,]),scale=F,center=T)))
#
#
##v <- vst(counts, method="anscombe.nb.simple", dispersion=dispersion)
#
##d <- NULL
#d <- design







