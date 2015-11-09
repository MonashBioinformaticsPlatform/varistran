
#
# All of the plots.
#

#' @export
shiny_filter <- function(y, counts=NULL, units="log2 counts", prefix="") {
    p <- function(name) paste0(prefix,name)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)

    ui <- shiny::tags$div(
        shiny::tags$h3("Select samples"),
        shiny::uiOutput(p("sample_selector")),
        shiny::tags$h3("Filter features"),
        shiny::numericInput(p("min_mean"), "Minimum mean expression level", 0.0),
        shiny::textOutput(p("report"))
    )

    server <- function(env) {
        env$output[[p("sample_selector")]] <- shiny::renderUI({
            y_val <- y(env)
            choices <- seq_len(ncol(y_val))
            if (!is.null(colnames(y_val)))
                names(choices) <- colnames(y_val)
            shiny::selectInput(p("samples"), "Select samples", selected=choices, choices=choices, multiple=TRUE)
        })

        env[[p("filtered")]] <- shiny::reactive({
            y_val <- y(env)
            counts_val <- counts(env)
            sample_select <- as.numeric(env$input[[p("samples")]])
            feature_select <- which(
                rowMeans(y_val[,sample_select,drop=FALSE]) >= env$input[[p("min_mean")]]
            )
            list(
                sample_select=sample_select,
                feature_select=feature_select,
                y=y_val[feature_select,sample_select,drop=FALSE],
                counts=counts_val[feature_select,sample_select,drop=FALSE]
            )
        })

        env$output[[p("report")]] <- shiny::renderText({
            filtered <- env[[p("filtered")]]()
            sprintf("%d of %d features will be used.",
                length(filtered$feature_select),
                nrow(y(env))
            )
        })
    }

    composable_shiny_app(ui, server)
}

#' Shiny report.
#'
#' Produce an interactive Shiny report showing diagnostic plots of transformed counts.
#'
#' @param y A matrix of exprssion levels, such as a transformed counts matrix.
#'
#' @param counts Optional. Original counts.
#'
#' @param prefix Optional, to fascilitate use as a component of a larger Shiny app. Inputs and outputs are given this prefix.
#'
#' @return A shiny.appobj.
#'
#' Used interactively, the shiny report runs immediately. Otherwise it can be launched by print()-ing it. A call to this function can also be the last line in an app.R file in a Shiny app directory.
#'
#' @author Paul Harrison
#'
#' @examples
#'
#' y <- varistran::vst(counts)
#' varistran::shiny_report(y, counts)
#'
#' @export
shiny_report <- function(y, counts=NULL, prefix="") {
    p <- function(name) paste0(prefix,name)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)

    filter <- shiny_filter(y, counts, prefix=p("filter_"))
    fy <- function(env) env[[p("filter_filtered")]]()$y
    fcounts <- function(env) env[[p("filter_filtered")]]()$counts

    stability <- shiny_stability(fy, fcounts, prefix=p("stability_"))
    biplot <- shiny_biplot(fy, prefix=p("biplot_"))
    heatmap <- shiny_heatmap(fy, prefix=p("heatmap_"))

    ui <- shiny::div(
        shiny::div(
            style="font-size: 150%; color: #bbbbbb; text-align: right; letter-spacing: 0.25em;",
            "Varistran"
        ),
        shiny::navlistPanel(
            widths=c(2,10),
            well=FALSE,
            shiny::tabPanel("Select and filter", filter$component_ui),
            shiny::tabPanel("Stability", stability$component_ui),
            shiny::tabPanel("Biplot", biplot$component_ui),
            shiny::tabPanel("Heatmap", heatmap$component_ui)
        )
    )

    server <- function(env) {
        filter$component_server(env)
        stability$component_server(env)
        biplot$component_server(env)
        heatmap$component_server(env)
    }

    composable_shiny_app(ui, server)
}
