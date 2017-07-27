
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
           varistran::vst(x, method="anscombe.nb.simple",
               dispersion=dispersion, cpm=T, design=design_for_estimation)
       ,
       "Anscombe" =
           varistran::vst(x, method="anscombe.nb",
               dispersion=dispersion, cpm=T, design=design_for_estimation)
       ,
       "Naive" =
           varistran::vst(x, method="naive.nb",
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
            name = rep(name, n),
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
    cat("\n\nArab\n\n")
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
    cat("\n\nBottomly\n\n")
    bottomly.eset <- load_bottomly()

    counts <- exprs(bottomly.eset)
    counts <- counts[rowSums(counts) >= 1, ]
    strain <- phenoData(bottomly.eset)$strain
    strain <- gsub("/", "_", strain)
    experiment.number <- factor( phenoData(bottomly.eset)$experiment.number )
    design <- model.matrix(~ strain + experiment.number)

    dds <- DESeqDataSetFromMatrix(
        counts,
        colData=data.frame(
            name=colnames(counts),
            strain=strain,
            experiment.number=experiment.number
        ),
        design = ~ strain + experiment.number)
    dds <- DESeq(dds)
    cat("DESeq2 dispersion function:\n")
    print(dds@dispersionFunction)

    real_stability_plot("stability_bottomly", counts, design, y_limit=0.3)
}

# Simulated data

if (T) {
    cat("\n\nSimulated\n\n")

    dispersion <- 0.1
    means <- 10 ^ seq(from=0,to=4,length.out=100000)
    nsamples <- 4

    # Always produce same output
    set.seed(2015)

    # Matrix update in place doesn't work, so use vector.
    # library(pryr)
    counts <- numeric(length(means) * nsamples)
    for(i in seq_along(means)) {
        counts[((i-1)*nsamples+1):(i*nsamples)] <-
            rnbinom(nsamples, size=1/dispersion, mu=means[i])

        # if (i%%1000 == 0) print(c(address(counts), refs(counts)))
    }

    counts <- matrix(counts, byrow=T,ncol=4)
    colnames(counts) <- c("a","b","c","d")

    design <- matrix(1, ncol=1,nrow=nsamples)

    dds <- DESeqDataSetFromMatrix(
        counts,
        colData=data.frame(names=as.character(seq_len(nsamples))),
        design=~1)
    dds <- DESeq(dds)
    cat("DESeq2 dispersion function:\n")
    print(dds@dispersionFunction)

    dgelist <- DGEList(counts)
    dgelist <- calcNormFactors(dgelist)
    dgelist <- estimateGLMCommonDisp(dgelist, design)

    dispersion <- dgelist$common.dispersion

    cat("EdgeR estimates the dispersion as",dispersion,"\n")


    stability_plot("stability_simulated", counts, design, dispersion, y_limit=0.55)
    stability_plot("stability_simulated_est", counts, design, NULL, y_limit=0.55)
    stability_plot("stability_simulated_est_blind", counts, design, NULL, y_limit=0.55, blind=T)
}
