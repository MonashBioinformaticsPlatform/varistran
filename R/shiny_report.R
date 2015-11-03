
#
# All of the plots.
#

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
shiny_report <- function(y, counts=NULL, prefix="") {
    p <- function(name) paste0(prefix,name)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)

    stability <- shiny_stability(y, counts, prefix=p("stability"))
    biplot <- shiny_biplot(y, prefix=p("biplot"))

    ui <- shiny::navlistPanel(
        widths=c(2,10),
        "Varistran",
        shiny::tabPanel("Stability", stability$component_ui),
        shiny::tabPanel("Biplot", biplot$component_ui)
    )

    server <- function(env) {
        stability$component_server(env)
        biplot$component_server(env)
    }

    composable_shiny_app(ui, server)
}
