
#
# Utility functions for creating composable shiny reports
#

#' @export
parenthetically <- function(...) {
    shiny::div(style="margin-top: 5em; font-size: 75%", ...)
}

#
# Shiny apps take most arguments as functions,
#   which may access reactive objects.
# Promote an argument to a function if necessary.
#

#' @export
ensure_reactable <- function(item) {
    if (is.function(item))
        item
    else
        function(env) item
}

#' Composable Shiny App.
#'
#' Create a Shiny app that can be used immediately, or incorporated into a larger Shiny app.
#'
#' Unlike a normal Shiny app, the server function is passed an environment, "env". This contains $input, $output and $server elements conventionally passed to a Shiny server function. However env can also be used to store reactive expressions that other components of a larger app might need to access.
#'
#' @param ui UI part of the Shiny app.
#'
#' @param server Server part of the Shiny app. This should take a single parameter, env.
#'
#' @return A shiny.appobj.
#'
#' This object has two additional properties, $component_ui and $component_server.
#'
#' A larger composable Shiny app would incorporate $component_ui into its own ui, and make sure to call $component_server(env) from from its own server function.
#'
#' @author Paul Harrison
#'
#' @export
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

#' Shiny plot.
#'
#' Convert a function that produces a plot into a Shiny app.
#'
#' @param callback Function to produce the plot. Takes one argument, "env" (see composable_shiny_app).
#'
#' @return A composable shiny.appobj.
#'
#' @author Paul Harrison
#'
#' @export
shiny_plot <- function(callback, width=500, height=500, dlname="plot", prefix="") {
    p <- function(name) paste0(prefix,name)

    ui <- shiny::tags$div(
        shiny::fluidRow(
          shiny::column(3, shiny::numericInput(p("width"), "Plot width", width, min=100, max=10000, step=50)),
          shiny::column(3, shiny::numericInput(p("height"), "Plot height", height, min=100, max=10000, step=50)),
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
