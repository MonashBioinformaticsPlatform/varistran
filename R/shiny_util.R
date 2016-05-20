
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

#' Make a value reactive if it is not already reactive.
#'
#' This is used to allow parameters to Shiny functions to be optionally reactive.
#'
#' @export
ensure_reactive <- function(item) {
    if ("reactive" %in% class(item))
        item
    else
        reactive(item)
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
#' Note that Shiny now has its own system for creating Shiny "modules". Varistran components will be converted to also support this style over time. shiny_plot already supports this.
#'
#' @author Paul Harrison
#'
#' @export
composable_shiny_app <- function(ui,server, title=NULL) {
    app <- shiny::shinyApp(
        shiny::fluidPage(ui, title=title),
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



# 
# Shiny plot
#

#' Shiny plot.
#'
#' Convert a function that produces a plot into a Shiny app.
#'
#' shiny_plot is the varistran-style composable app object. shiny_plot_ui and shiny_plot_server are appropriate for the normal Shiny module system.
#'
#' @param callback Function to produce the plot. In shiny_plot, this takes one argument, "env" (see composable_shiny_app). In shiny_plot_server, it takes no arguements.
#'
#' @param ... Extra arguments to shiny_plot and shiny_plot_ui are passed on to plotOutput, for example brushing options.
#'
#' @return shiny_plot() returns a composable shiny.appobj.
#'
#' @author Paul Harrison
#'
#' @export
shiny_plot <- function(callback, width=500, height=500, dlname="plot", prefix="", ...) {
    composable_shiny_app(
        shiny_plot_ui(prefix, width=width, height=height, ...),
        function(env) {
            callback_no_env <- function() callback(env) 
            callModule(shiny_plot_server, prefix, callback_no_env, dlname, session=env$session)
        }
    )
}


#' @describeIn shiny_plot UI part of Shiny module.
#' @export
shiny_plot_ui <- function(id, width=500, height=500, ...) {
    ns <- NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
          shiny::column(3, shiny::numericInput(ns("width"), "Plot width", width, min=100, max=10000, step=50)),
          shiny::column(3, shiny::numericInput(ns("height"), "Plot height", height, min=100, max=10000, step=50)),
          shiny::column(4, shiny::tags$label("Download"), shiny::tags$br(),
              shiny::downloadButton(ns("pdf"), "PDF"),
              shiny::downloadButton(ns("eps"), "EPS"),
              shiny::downloadButton(ns("png"), "PNG"))
        ),
        shiny::plotOutput(ns("plot"), width="auto", height="auto", ...)
    )
}

#' @describeIn shiny_plot server part of Shiny module, to be invoked from a server function with callModule.
#' @export
shiny_plot_server <- function(input, output, session, callback, dlname="plot") {
    output$plot <- shiny::renderPlot(
        { callback() },
        width=function() input$width,
        height=function() input$height
    )
    
    output$pdf <- shiny::downloadHandler(
        paste0(dlname,".pdf"),
        function(filename) {
            pdf(filename, width=input$width/72, height=input$height/72)
            callback()
            dev.off()
        }
    )
    
    output$eps <- shiny::downloadHandler(
        paste0(dlname,".eps"),
        function(filename) {
            postscript(filename, width=input$width/72, height=input$height/72,
                       paper="special", onefile=FALSE, horizontal=FALSE)
            callback()
            dev.off()
        }
    )
    
    output$png <- shiny::downloadHandler(
        paste0(dlname,".png"),
        function(filename) {
            png(filename, width=input$width, height=input$height)
            callback()
            dev.off()
        }
    )
}

