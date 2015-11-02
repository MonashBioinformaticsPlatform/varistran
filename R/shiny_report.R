
#
# All of the plots.
#

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
