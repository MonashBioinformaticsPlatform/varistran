
#
# All of the plots.
#

#' @export
shiny_filter <- function(y, counts=NULL, units="log2 counts", prefix="") {
    p <- function(name) paste0(prefix,name)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)

    ui <- shiny::tags$div(
        titlePanel("Select and filter"),
        numericInput(p("min_mean"), "Minimum mean expression level", 0.0),
        textOutput(p("report"))
    )

    server <- function(env) {
        env[[p("filtered")]] <- reactive({
            y_val <- y(env)
            counts_val <- counts(env)
            select <- rowMeans(y_val) >= env$input[[p("min_mean")]]
            list(
                select=select,
                y=y_val[select,,drop=FALSE],
                counts=counts_val[select,,drop=FALSE]
            )
        })

        env$output[[p("report")]] <- renderText({
            filtered <- env[[p("filtered")]]()
            sprintf("%d of %d features will be used.",
                sum(filtered$select),
                length(filtered$select)
            )
        })
    }

    composable_shiny_app(ui, server)
}

#' Shiny report.
#'
#' Produce an interactive Shiny report showing diagnostic plots of transformed counts.
#'
#' @param y Transformed counts.
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
shiny_report <- function(y, counts=NULL, units="log2 count", prefix="") {
    p <- function(name) paste0(prefix,name)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)

    filter <- shiny_filter(y, counts, prefix=p("filter_"))
    fy <- function(env) env[[p("filter_filtered")]]()$y
    fcounts <- function(env) env[[p("filter_filtered")]]()$counts

    stability <- shiny_stability(fy, fcounts, prefix=p("stability_"))
    biplot <- shiny_biplot(fy, prefix=p("biplot_"))
    heatmap <- shiny_heatmap(fy, units=units, prefix=p("heatmap_"))

    ui <- shiny::navlistPanel(
        widths=c(2,10),
        "Varistran",
        shiny::tabPanel("Select and filter", filter$component_ui),
        shiny::tabPanel("Stability", stability$component_ui),
        shiny::tabPanel("Biplot", biplot$component_ui),
        shiny::tabPanel("Heatmap", heatmap$component_ui)
    )

    server <- function(env) {
        filter$component_server(env)
        stability$component_server(env)
        biplot$component_server(env)
        heatmap$component_server(env)
    }

    composable_shiny_app(ui, server)
}
