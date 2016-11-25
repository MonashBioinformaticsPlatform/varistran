
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
ensure_reactable <- function(item, env_name=NULL) {
    if (!is.null(env_name) && is.null(item))
        function(env) env[[env_name]]()
    else if (is.function(item))
        item
    else
        function(env) item
}


#' Make a value reactive if it is not already reactive.
#'
#' This is used to allow parameters to Shiny functions to be optionally reactive.
#'
#' @export
ensure_reactive <- function(item, name=NULL, env=NULL, default=function() stop("Missing ",name)) {
    if (is.null(item) && !is.null(env))
        shiny::reactive({ 
            if (is.null(env[[name]]))
                default()
            else
                env[[name]]() 
        })
    else if ("reactive" %in% class(item))
        item
    else
        shiny::reactive(item)
}


#' Call object with "request" if it is callable.
#'
#' This is used to support older UI code which doesn't wrap UI in function(request) { ... }.
#'
#' @export
call_ui <- function(ui, request) {
    if (!is.function(ui)) {
        if (!is.character(ui))
            warning("Using UI not wrapped in function. Bookmarks may not restore.")
        return(ui)
    }
    
    ui(request)
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
#' @param enableBookmarking Parameter passed to shinyApp. A bookmark button is displayed unless enableBookmarking is "disable".
#'
#' @param output Name in env that will contain a reactive value that produces an output (eg for use as a gadget).
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
composable_shiny_app <- function(ui, server, title=NULL, enableBookmarking="server", output=NULL) {
    app <- shiny::shinyApp(
        function(request) shiny::fluidPage(
            call_ui(ui, request), 
            if (enableBookmarking != "disable") 
                shiny::div(
                    style="position: fixed; bottom: 5px; left: 5px;", 
                    shiny::bookmarkButton("bookmark")),
            shiny::div(style="height: 5em"),
            title=title),
        function(input,output,session) {
            # Create a shared environment to
            # stash reactive values in.
            env <- new.env()
            env$input <- input
            env$output <- output
            env$session <- session
            server(env)
        },
        
        enableBookmarking = enableBookmarking
    )
    app$component_title <- title
    app$component_ui <- ui
    app$component_server <- server
    app$component_output <- output
    app
}



#' Run a composable shiny app as a gadget
#'
#' @export
run_as_gadget <- function(app, viewer=shiny::dialogViewer("Gadget")) {
    shiny::runGadget(
        function(request) miniUI::miniPage(
            miniUI::gadgetTitleBar(app$component_title),
            call_ui(app$component_ui, request)),
        
        function(input, output, session) {
            env <- new.env()
            env$input <- input
            env$output <- output
            env$session <- session
            app$component_server(env)
            
            shiny::observeEvent(env$input$done, {
                result <- 
                    if (!is.null(app$component_output))
                        env[[app$component_output]]()
                    else
                        NULL
                
                shiny::stopApp(result)
            })
        },
        
        viewer=viewer    
    )
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
shiny_plot <- function(callback=NULL, width=500, height=500, dlname="plot", prefix="", ...) {
    composable_shiny_app(
        shiny_plot_ui(prefix, width=width, height=height, ...),
        function(env) {
            if (!is.null(callback))
                env[[paste0(prefix,"-callback")]] <- function() callback(env)
            callback_no_env <- function() env[[paste0(prefix,"-callback")]]()
            shiny::callModule(shiny_plot_server, prefix, callback_no_env, dlname, session=env$session)
        }
    )
}


#' @describeIn shiny_plot UI part of Shiny module.
#' @export
shiny_plot_ui <- function(id, width=500, height=500, ...) {
    ns <- shiny::NS(id)

    function(request) shiny::tagList(
        shiny::fluidRow(
          shiny::column(3, shiny::numericInput(ns("width"), "Plot width", width, min=100, max=10000, step=50)),
          shiny::column(3, shiny::numericInput(ns("height"), "Plot height", height, min=100, max=10000, step=50)),
          shiny::column(6, shiny::tags$label("Download"), shiny::tags$br(),
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
