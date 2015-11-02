
#
# Utility functions for creating composable shiny reports
#

#
# Shiny apps take most arguments as functions,
#   which may access reactive objects.
# Promote an argument to a function if necessary.
#
ensure_reactable <- function(item) {
    if (is.function(item))
        item
    else
        function(env) item
}


composable_shiny_app <- function(ui,server) {
    app <- shiny::shinyApp(
        shiny::fluidPage(ui),
        function(input,output,session) {
            # Create a shared environment to
            # stash reactive values in.
            env <- new.env()
            env$input <- input
            env$output <- output
            env$session <- session
            server(env)
        }
    )
    app$component_ui <- ui
    app$component_server <- server
    app
}


shiny_plot <- function(callback, width=500, height=500, dlname="plot", prefix="") {
    p <- function(name) paste0(prefix,name)

    ui <- shiny::tags$div(
        shiny::fluidRow(
          shiny::column(3, shiny::numericInput(p("width"), "Plot width", width, min=100, max=10000, step=50)),
          shiny::column(3, shiny::numericInput(p("height"), "Plot height", height, width, min=100, max=10000, step=50)),
          shiny::column(4, shiny::tags$label("Download"), shiny::tags$br(),
              shiny::downloadButton(p("pdf"), "PDF"),
              shiny::downloadButton(p("eps"), "EPS"))
        ),
        shiny::plotOutput(p("plot"), width="auto", height="auto")
    )

    server <- function(env) {
        output <- env$output
        i <- function(name) env$input[[p(name)]]

        output[[p("plot")]] <- shiny::renderPlot(
            { callback(env) },
            width=function() i("width"),
            height=function() i("height")
        )

        output[[p("pdf")]] <- shiny::downloadHandler(
            paste0(dlname,".pdf"),
            function(filename) {
                pdf(filename, width=i("width")/72, height=i("height")/72)
                callback(env)
                dev.off()
            }
        )

        output[[p("eps")]] <- shiny::downloadHandler(
            paste0(dlname,".eps"),
            function(filename) {
                postscript(filename, width=i("width")/72, height=i("height")/72,
                           paper="special", onefile=FALSE, horizontal=FALSE)
                callback(env)
                dev.off()
            }
        )
    }

    composable_shiny_app(ui, server)
}
