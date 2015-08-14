
source("test/common.R")

library("MASS")
library("edgeR")
library("DESeq2")


window <- function(vec, size) {
    cs <- c(0.0, cumsum(vec))
    (cs[(1+size):length(cs)] - cs[1:(length(cs)-size)]) / size
}

bin <- function(vec, n) {
    result <- numeric(n)
    for(i in seq_len(n))
        result[i] <- mean(vec[ (floor((i-1)*length(vec)/n)+1):floor(i*length(vec)/n) ])
    result
}

#library("NBPSeq")
#data("arab")
#good <- rowSums(arab) > 0
#counts <- arab[good,]
#group <- factor(c("a","a","a","b","b","b"))
#batch <- factor(c("b1","b2","b3","b1","b2","b3"))
#design <- model.matrix(~ group + batch)


#dat <- read.csv("/data/home/pfh/2014/gld2combined/pipeline/raw/genewise-count.csv",row.names=1)
#dat <- dat[,1:9]
#good <- rowSums(dat) > 0
#counts <- dat[good,]
#group <- factor(c("N2","N2","N2","Gld2","Gld2","Gld2","Cpb3","Cpb3","Cpb3"),
#    levels=c("N2","Gld2","Cpb3"))
#design <- model.matrix(~ group)


dispersion <- 0.01
means <- 2 ^ seq(from=0,to=16,length.out=10000)
counts <- matrix(0.0, nrow=length(means), ncol=4)
colnames(counts) <- c("a","b","c","d")
group <- factor(c("1","1","2","2"))
for(i in seq_len(length(means))) {
    counts[i,] <- rnbinom(ncol(counts), size=1/dispersion, mu=means[i])
}
design.formula <- ~ group
design <- model.matrix(design.formula)

tnull <- t(Null(design))
y <- DGEList(counts)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)

means <- rowMeans(counts)
totals <- rowSums(counts)

dispersion <- y$common.dispersion
#dispersion <- min(y$trended.dispersion)
print(dispersion)


#dispersion <- 1.0 #dispersion * 0.4
print(dispersion)

#dispersion <- optimal.dispersion(counts, design=design)
#print(dispersion)

#dispersion <- NULL

#dds <- DESeqDataSetFromMatrix(counts, colData=data.frame(name=colnames(counts),group=group), #design = design.formula)
#dds <- DESeq(dds)
#print(dds@dispersionFunction)


#v <- vst(counts,method="anscombe.nb",dispersion=dispersion,cpm=T)
#plotMDS(v, top=nrow(v))
#
#low <- apply(v,1,min)
#high <- apply(v,1,max)
#span <- high-low
#good <- order(span, decreasing=T)[1:500]
#nesoni.heatmap(t(scale(t(v[good,]),scale=F,center=T)))


#v <- vst(counts, method="anscombe.nb.simple", dispersion=dispersion)

#d <- NULL
d <- design

vs = list(
   "Naive" = vst(counts, method="naive.nb", dispersion=dispersion, cpm=T, design=d)
   ,
   "log(x+c)" = vst(counts, method="anscombe.nb.simple", dispersion=dispersion, cpm=T, design=d)
   ,
   "Anscombe" = vst(counts, method="anscombe.nb", dispersion=dispersion, cpm=T, design=d)
   
   #cpm = cpm(y, prior.count=0.5/dispersion, log=T)
)

ord <- order(totals)

dfs <- list()

for(name in names(vs)) {
    print(name)
    v <- vs[[name]]
    print(.dispersion.score(v, design))
    variances <- apply(v,1,function(x) {
       r <- tnull %*% x
       mean(r*r)
    })

    
    df <- data.frame(
        name=rep(name, length(ord)),
        rank=seq_len(length(ord)),
        means=means[ord],
        variances = variances[ord],
        sds = sqrt(variances[ord]))
    
    #dfs[[length(dfs)+1]] <- df
    
    #r <- 100
    #windowed <- df[r:(nrow(df)-r),]
    #for(i in 1:nrow(windowed))
    #    windowed[i,"variances"] <- mean(df[i:(i+r*2+1),"variances"])
    
    #n <- 1000
    #windowed <- data.frame(
    #    rank = window(df$rank, n),
    #    totals = window(df$totals, n),
    #    variances = window(df$variances, n)
    #    )
    #windowed$name <- rep(name,nrow(windowed))
    ##windowed$variances <- windowed$variances / mean(windowed$variances)
    #    
    #dfs[[length(dfs)+1]] <- windowed


    n <- 20
    binned <- data.frame(
        rank = bin(df$rank, n),
        means = bin(df$means, n),
        variances = bin(df$variances, n)
        )
    
    binned$name <- rep(name,nrow(binned))
    dfs[[length(dfs)+1]] <- binned
}

bigdf <- do.call(rbind, dfs)

to <- floor(log10(max(bigdf$means)))

print(
  ggplot(bigdf, aes(x=means,y=sqrt(variances),color=name)) 
  + geom_line() 
  + geom_point(size=3) 
  + scale_x_log10(breaks=10^(0:to))
  + ylim(0.0, NA)
  + labs(x="Mean count", y="sqrt(variance of transformed counts)", color="Transformation")
  + theme_bw()
  + theme(legend.justification=c(1,0), legend.position=c(1,0))
  )
  
#plot(rank(totals,ties.method="random"),sqrt(variances),cex=0.1)






