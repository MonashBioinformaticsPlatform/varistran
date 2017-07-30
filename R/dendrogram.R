
# This code for drawing dendrograms has been adapted from ggdendro package.

dendro_data <- function(model, type=c("rectangle","triangle"), ...) {
  dhc <- as.dendrogram(model)
  hcdata <- dendrogram_data(dhc, type=type, ...)
  list(
      segments = hcdata$segments,
      labels = hcdata$labels,
      class="hclust"
  )
}

dendrogram_data <- function (x, type = c("rectangle", "triangle"), ...){ 

    # Initialise variables that used to be in parameter list
    leaflab <- "perpendicular"
    center <- FALSE
    xlab <- ""
    ylab <- ""
    horiz <- FALSE
    #frame.plot <- FALSE
    xaxt <- "n"
    yaxt <- "s"
    dLeaf <- NULL 
    edge.root <- is.leaf(x) || !is.null(attr(x, "edgetext"))

    type <- match.arg(type)
    #leaflab <- match.arg(leaflab)
    hgt <- attr(x, "height")
    if (edge.root && is.logical(edge.root)) 
        edge.root <- 0.0625 * if (is.leaf(x)) 1 else hgt
    mem.x <- .memberDend(x)
    yTop <- hgt + edge.root
    if (center) {
        x1 <- 0.5
        x2 <- mem.x + 0.5
    }
    else {
        x1 <- 1
        x2 <- mem.x
    }
    xl. <- c(x1 - 1/2, x2 + 1/2)
    yl. <- c(0, yTop)
    if (edge.root) {
        if (!is.null(et <- attr(x, "edgetext"))) {
            my <- mean(hgt, yTop)
        }
    }
    
    gg.plotNode <- function (x1, x2, subtree, type, center, leaflab, dLeaf, horiz=FALSE, ddsegments=NULL, ddlabels=NULL) {
        inner <- !is.leaf(subtree) && x1 != x2
        yTop <- attr(subtree, "height")
        bx <- plotNodeLimit(x1, x2, subtree, center)
        xTop <- bx$x
        hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
        Xtract <- function(nam, L, default, indx) rep(if (nam %in% 
                                    names(L)) L[[nam]] else default, length.out = indx)[indx]
        asTxt <- function(x){
            if (is.character(x) || is.expression(x))
                x else
            if (is.null(x)) "" else as.character(x)
        }
        i <- if (inner || hasP) 
                    1
                else 2
        lab.cex <- 1
        if (is.leaf(subtree)) {
            if (leaflab == "perpendicular") {
                Y <- yTop - dLeaf * lab.cex
                X <- xTop
                srt <- 90
                adj <- 1
                nodeText <- asTxt(attr(subtree, "label"))
                ddlabels <- rbind(ddlabels, data.frame(x=X, y=0, text=nodeText))
            }
        }
        else if (inner) {
            segmentsHV <- function(x0, y0, x1, y1) {
                # *************************
                data.frame(x0, y0, x1, y1) #AdV
            }
            for (k in seq_along(subtree)) {
                child <- subtree[[k]]
                yBot <- attr(child, "height")
                if (getOption("verbose")) 
                    cat("ch.", k, "@ h=", yBot, "; ")
                if (is.null(yBot)) 
                    yBot <- 0
                xBot <- if (center) 
                            mean(bx$limit[k:(k + 1)])
                        else bx$limit[k] + .midDend(child)
                if (type == "triangle") {
                    ddsegments <- rbind(ddsegments, segmentsHV(xTop, yTop, xBot, yBot))
                }
                else {
                    ddsegments <- rbind(ddsegments, segmentsHV(xTop, yTop, xBot, yTop))
                    ddsegments <- rbind(ddsegments, segmentsHV(xBot, yTop, xBot, yBot))
                }
                vln <- NULL
                if (!is.null(attr(child, "edgetext"))) {
                    edgeText <- asTxt(attr(child, "edgetext"))
                    if (!is.null(vln)) {
                        mx <- if (type == "triangle") 
                                    (xTop + xBot + ((xTop - xBot)/(yTop - yBot)) * vln)/2
                                else xBot
                        my <- (yTop + yBot + 2 * vln)/2
                    }
                    else {
                        mx <- if (type == "triangle") 
                                    (xTop + xBot)/2
                                else xBot
                        my <- (yTop + yBot)/2
                    }
                }
                plotNode_result <- gg.plotNode(bx$limit[k], bx$limit[k + 1], subtree = child, 
                        type, center, leaflab, dLeaf, horiz, ddsegments, ddlabels)
                ddsegments <- plotNode_result$segments
                ddlabels <- plotNode_result$labels
            }
        }
        return(list(segments=ddsegments, labels=ddlabels))
    }
    
    ret <- gg.plotNode(x1, x2, x, type = type, center = center, leaflab = leaflab, 
            dLeaf = dLeaf, horiz=FALSE, 
            ddsegments=NULL, ddlabels=NULL)
    names(ret$segments) <- c("x", "y", "xend", "yend")
    names(ret$labels) <- c("x", "y", "label")
    ret
}



.memberDend <- function (x) 
{
  r <- attr(x, "x.member")
  if (is.null(r)) {
    r <- attr(x, "members")
    if (is.null(r)) 
      r <- 1L
  }
  r
}


plotNodeLimit <- function (x1, x2, subtree, center) 
{
  inner <- !is.leaf(subtree) && x1 != x2
  if (inner) {
    K <- length(subtree)
    mTop <- .memberDend(subtree)
    limit <- integer(K)
    xx1 <- x1
    for (k in 1L:K) {
      m <- .memberDend(subtree[[k]])
      xx1 <- xx1 + (if (center) 
        (x2 - x1) * m/mTop
                    else m)
      limit[k] <- xx1
    }
    limit <- c(x1, limit)
  }
  else {
    limit <- c(x1, x2)
  }
  mid <- attr(subtree, "midpoint")
  center <- center || (inner && !is.numeric(mid))
  x <- if (center) 
    mean(c(x1, x2))
  else x1 + (if (inner) 
    mid
             else 0)
  list(x = x, limit = limit)
}


.midDend <- function (x) 
  if (is.null(mp <- attr(x, "midpoint"))) 0 else mp



