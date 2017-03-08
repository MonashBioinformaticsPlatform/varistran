
#
# grid graphics utility functions
#

# This is necessary because we are providing new methods for S3 generics
# defined by grid.
#
#' @import grid


# I've had problems if varistran functions are invoked directly with 
# varistran::... without library being called first.
#
# However, the problem no long seems to be occurring.
#
#.onLoad <- function(...) {
#    library(grid)
#}


.default <- function(a,b) if (is.null(a)) b else a

#
# varistran grid class
# - print()-able
# - can have manually specified size, for layout purposes
#

#' @export
heightDetails.varistran_grob <- function(x) {
    .default(x$height_hint, grid::unit(1,"null"))
}

#' @export
widthDetails.varistran_grob <- function(x) {
    .default(x$width_hint, grid::unit(1,"null"))
}

#' @export
print.varistran_grob <- function(x, newpage=TRUE, ...) {
    if (newpage)
        grid::grid.newpage()
    grid::grid.draw(x)
}

#' Make a grob better.
#'
#' Give a grob a margin, and size hints for layout, and make it print()-able.
#'
#' @param grob to wrap.
#'
#' @param width Size hint.
#' 
#' @param height Size hint.
#'
#' @param pad Padding to place around \code{grob}.
#'
#' @export
varistran_grob <- function(
        grob,
        width=unit(1,"null"),
        height=unit(1,"null"),
        pad=0) {
    if (is.numeric(pad))
        pad <- unit(pad, "lines")

    # Null grobs take no space
    if ("null" %in% class(grob))
        pad <- unit(0,"lines")

    if (identical(height,"inherit"))
        height <- grobHeight(grob) + 2*pad
    if (identical(width,"inherit"))
        width <- grobWidth(grob) + 2*pad

    gTree(
        children=gList(grob),
        cl="varistran_grob",
        width_hint=width,
        height_hint=height,
        vp=viewport(
            width=unit(1,"npc")-2*pad,
            height=unit(1,"npc")-2*pad
        )
    )
}

#' Grob to display an ordering as a dendrogram.
#'
#' Display the dendrogram of an ordering, as a grid grob.
#'
#' @param ordering An ordering produced by the \code{make_ordering} function.
#'
#' @param transpose Show a horizontal rather than a vertical dendrogram.
#'
#' @param mirror Should the leaves point left rather than right (down rather than up if \code{transpose=TRUE})?
#'
#' @param hint_size Width (height if \code{transpose=TRUE}).
#'
#' @seealso \code{\link{make_ordering}}
#'
#' @export
ordering_grob <- function(ordering, transpose=FALSE, mirror=FALSE, hint_size=unit(5,"lines")) {
    if (is.null(ordering$dendrogram))
        return(nullGrob())

    dd <- ggdendro::dendro_data(ordering$dendrogram)
    dds <- dd$segments
    x0 <- dd$segments$x
    x1 <- dd$segments$xend
    y0 <- dd$segments$y
    y1 <- dd$segments$yend
    xscale <- range(x0,x1) + c(-0.5,0.5)
    yscale <- range(y0,y1)
    width <- unit(1,"null")
    height <- hint_size

    # Allow zero width dendrogram
    if (yscale[1] == yscale[2])
        yscale <- yscale + c(-1,1)

    if (mirror) {
        yscale <- c(yscale[2],yscale[1])
    }

    if (transpose) {
        t0 <- x0
        t1 <- x1
        tscale <- xscale
        tsize <- width
        x0 <- y0
        x1 <- y1
        xscale <- yscale
        width <- height
        y0 <- t0
        y1 <- t1
        yscale <- tscale
        height <- tsize
    }

    segments <- segmentsGrob(
        x0=x0, y0=y0, x1=x1, y1=y1, default.units="native",
        vp=viewport(xscale=xscale, yscale=yscale),
        )

    varistran_grob(segments, width=width, height=height)
}


#' A palette from cool to hot.
#'
#' @export
unsigned_colors <- hsv(
    h=seq(0.95,1.15, length.out=256)%%1.0,
    v=seq(0,1, length.out=256)**0.5,
    s=seq(1,0,length.out=256)**0.5)


#' A red-blue palette suitable for signed data.
#'
#' @export
signed_colors <- hsv(
    h=(sign(seq(-1.0,1.0, length.out=256))*0.2+0.8)%%1.0,
    v=1,
    s=abs(seq(-1,1,length.out=256)))

#' A grob for heatmaps.
#'
#' @param data Matrix of data to display.
#'
#' @param signed Should a signed palette be used?
#'
#' @param legend_title Title for legend grob.
#'
#' @param vp_name Viewport name for heatmap grob. This is useful if you will need to refer to the grob later, for example for interaction.
#'
#' @return A list containing a heatmap grob and associated legend grob.
#'
#' @export
heatmap_grob <- function(data, signed=TRUE, legend_title="", vp_name=NULL) {
    if (signed) {
        radius <- max(abs(data), na.rm=TRUE)
        range <- c(-radius, radius)
        col <- signed_colors
    } else {
        range <- c(0, max(data, na.rm=TRUE))
        col <- unsigned_colors
    }

    # Avoid divide by zero
    if (range[1] == range[2])
        range <- range + c(-1,1)

    scaled <- pmin(pmax(as.integer( (data-range[1])/(range[2]-range[1])*256+1),1),256)

    fill <- col[scaled]
    fill[ is.na(fill) ] <- "#888888"

    heatmap <- rectGrob(
        x=rep(seq_len(ncol(data))-1, each=nrow(data)),
        y=rep(seq_len(nrow(data))-1, ncol(data)),
        width=1,
        height=1,
        just=c(0,0),
        default.units="native",
        gp=gpar(col=NA, fill=fill),
        vp=viewport(xscale=c(0,ncol(data)),yscale=c(0,nrow(data)), name=vp_name)
    )

    legend_heatmap <- rectGrob(
        x=seq_along(col)-1,
        y=0,
        width=1,
        height=1,
        just=c(0,0),
        default.units="native",
        gp=gpar(col=NA,fill=col),
        vp=viewport(xscale=c(0,length(col)),yscale=c(0,1))
    )

    axis <- xaxisGrob(
            at=axisTicks(range,log=FALSE,nint=5),
            vp=viewport(width=1,height=0,y=7/8,xscale=range),
            gp=gpar(cex=0.75)
            )

    text <- textGrob(legend_title)

    legend <- frameGrob()
    legend <- packGrob(legend, varistran_grob(legend_heatmap,width=unit(5,"lines"),height=unit(1,"lines")), row=1)
    legend <- packGrob(legend, varistran_grob(axis,height=unit(2,"lines")), row=2)
    legend <- packGrob(legend, varistran_grob(text,height="inherit"), row=3)

    list(
        heatmap=heatmap,
        legend=legend
    )
}

#' A collection of horizontal labels which shrink if they get too small.
#'
#' @param label Vector of labels.
#'
#' @param x Vector of x positions.
#'
#' @param y Vector of y positions.
#'
#' @param just Text justification.
#'
#' @param max_width Maximum width of labels (defaults to units of lines).
#'
#' @param ... Further arguments to \code{grob}.
#'
#' @export
shrinktext_grob <- function(label,x,y,just,max_width=20,...) {
    if (is.numeric(max_width))
        max_width <- unit(max_width, "lines")

    grob(label=label,x=x,y=y,just=just,shrinktext_max_width=max_width,...,cl="shrinktext_grob")
}


#' @export
widthDetails.shrinktext_grob <- function(x) {
    min(x$shrinktext_max_width, max(stringWidth(x$label)))
}


#' @export
drawDetails.shrinktext_grob <- function(x, recording) {
    ratio <- convertHeight(unit(1,"native"), "lines", valueOnly=TRUE)
    cex <- min(1.0, ratio)

    if (cex < 0.2) return()

    grid.text(
        x$label,
        x=x$x,
        y=x$y,
        just=x$just,
        default.units="native",
        gp=gpar(cex=cex)
    )
}


#' A collection of vertical labels which shrink if they get too small.
#'
#' @param label Vector of labels.
#'
#' @param x Vector of x positions.
#'
#' @param y Vector of y positions.
#'
#' @param just Text justification.
#'
#' @param ... Further arguments to \code{grob}.
#'
#' @export
vertical_shrinktext_grob <- function(label,x,y,just,...) {
    grob(label=label,x=x,y=y,just=just,...,cl="vertical_shrinktext_grob")
}


#' @export
heightDetails.vertical_shrinktext_grob <- function(x) {
    max(stringWidth(x$label))
}


#' @export
drawDetails.vertical_shrinktext_grob <- function(x, recording) {
    ratio <- convertWidth(unit(1,"native"), "lines", valueOnly=TRUE)
    cex <- min(1.0, ratio)

    if (cex < 0.2) return()

    grid.text(
        x$label,
        x=x$x,
        y=x$y,
        just=x$just,
        default.units="native",
        rot=90,
        gp=gpar(cex=cex)
    )
}
