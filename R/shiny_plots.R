#
# Shiny wrappers around plotting functions
#

#' @export
shiny_stability <- function(y, x=NULL, design=NULL, bins=20, prefix="") {
    ns <- shiny::NS(prefix)

    y <- ensure_reactable(y)
    x <- ensure_reactable(x)
    design <- ensure_reactable(design)
    bins <- ensure_reactable(bins)

    plot <- shiny_plot(
        function(env) {
            bins <- env$input[[ns("bins")]]
            stopifnot(bins >= 1)
            print(plot_stability(
                y=y(env),
                x=x(env),
                design=design(env),
                bins=bins
            ))
        },
        dlname="stability_plot",
        prefix=ns("plot")
    )

    ui <- function(request) shiny::tags$div(
        shiny::h3("Stability plot"),
        shiny::numericInput(ns("bins"), "Bins", value=20,min=1,max=10000),
        call_ui(plot$component_ui, request),
        shiny::p(
            "This is a diagnostic plot for the variance stabilizing transformation. ",
            "It shows the amount of noise in the transformed data ",
            "across different expression levels. ",
            "If the transformation has been successful, ",
            "this should be close to a flat line."
            ),
        parenthetically("This plot is produced by varistran::plot_stability.")
    )

    server <- function(env) {
        plot$component_server(env)
    }

    composable_shiny_app(ui, server)
}


#' Shiny wrapper for limma's MDS plot
#'
#' @export
shiny_mds_plot <- function(x, sample_labels=NULL, prefix="") {
    ns <- shiny::NS(prefix)

    x <- ensure_reactable(x)
    sample_labels <- ensure_reactable(sample_labels)
    
    plot <- shiny_plot(
        function(env) {
            x_val <- x(env)
            sample_labels_val <- sample_labels(env)
            
            if (!is.null(sample_labels_val))
                colnames(x_val) <- sample_labels_val
        
            limma::plotMDS(x_val, env$input[[ns("genes")]])
        },
        prefix = paste0(prefix,"plot")
    )
    
    ui <- function(request) shiny::div(
        shiny::h3("limma MDS plot"),
        shiny::numericInput(ns("genes"), "Use this many top genes", value=500, min=1,max=20000),
        call_ui(plot$component_ui, request),
        parenthetically("This plot is produced by limma::plotMDS. Gene selection is \"pairwise\".")
    )
    
    server <- function(env) {
        plot$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}


#' @export
shiny_biplot <- function(x, sample_labels=NULL, feature_labels=NULL, n_features=20, balance=0.25, text_size=0.025, prefix="") {
    ns <- shiny::NS(prefix)

    x <- ensure_reactable(x)
    sample_labels <- ensure_reactable(sample_labels)
    feature_labels <- ensure_reactable(feature_labels)

    plot <- shiny_plot(
        function(env) {
            print(plot_biplot(
                x = x(env),
                sample_labels = sample_labels(env),
                feature_labels = feature_labels(env),
                n_features = env$input[[ns("n_features")]],
                balance = env$input[[ns("balance")]],
                text_size = 0.025
            ))
        },
        dlname="biplot",
        prefix=ns("plot")
    )

    ui <- function(request) shiny::tags$div(
        shiny::h3("Biplot"),
        shiny::numericInput(ns("n_features"), "Number of labelled features", 20, min=0, step=1),
        shiny::numericInput(ns("balance"), "Feature/sample relative scaling", 0.25, min=0,step=0.05),
        call_ui(plot$component_ui, request),
        parenthetically("This plot is produced by varistran::plot_biplot.")
    )

    server <- function(env) {
        plot$component_server(env)
    }

    composable_shiny_app(ui, server)
}

#' @export
shiny_heatmap <- function(y, sample_labels=NULL, feature_labels=NULL, prefix="") {
    ns <- shiny::NS(prefix)

    y <- ensure_reactable(y)
    sample_labels <- ensure_reactable(sample_labels)
    feature_labels <- ensure_reactable(feature_labels)

    plot <- shiny_plot(
        callback = function(env) shiny::withProgress(message="Plotting", {
            print(env[[ns("grob")]]())
            
            # Configure base graphics to match heatmap
            seekViewport("heatmap")
            pltvec <- gridBase::gridPLT()
            # Selection to span full width of plot
            pltvec <- c(0, 1, pltvec[3], pltvec[4])
            par(new=T, plt=pltvec)
            plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1),
                 xaxs="i", yaxs="i")
        }),
        width=700,
        height=700,
        dlname="heatmap",
        prefix=ns("plot"),
        brush=shiny::brushOpts(
            id=ns("brush"),  
            direction = 'y',
            resetOnNew = TRUE,
            clip = TRUE,
            delay = 600000
        )
    )

    ui <- function(request) shiny::tags$div(
        shiny::h3("Heatmap"),
        shiny::p("Features are selected based on span of expression levels."),
        shiny::numericInput(ns("n"), "Number of features to show", 50, min=10,max=2000,step=10),
        shiny::checkboxInput(ns("cluster_samples"), "Cluster samples", FALSE),
        call_ui(plot$component_ui, request),
        #shiny::uiOutput(ns("selected_text")),
        parenthetically("This plot is produced by varistran::plot_heatmap.")
    )

    server <- function(env) {
        e <- function(name) env[[ns(name)]]()
        
        y_val <- shiny::reactive( as.matrix(y(env)) )
    
        env[[ns("selection")]] <- shiny::reactive({
            n <- env$input[[ns("n")]]
            if (n > 2000) stop("Drawing large heatmaps uses excessive system resources. Sorry.")

            y_span <- apply(y_val(),1,max) - apply(y_val(),1,min)
            selection <- rep(FALSE,nrow(y_val()))
            selection[ order(-y_span)[ seq_len(n) ] ] <- TRUE
            selection <- which(selection)
            
            if (length(selection) < 1) stop("No features to show.")

            selection
        })
    
        env[[ns("grob")]] <- shiny::reactive({
            plot_heatmap(
                y=y_val()[e("selection"),,drop=FALSE],
                sample_labels=sample_labels(env),
                feature_labels=feature_labels(env)[e("selection")],
                cluster_samples=env$input[[ns("cluster_samples")]]
            )
        })

        plot$component_server(env)
        
        env[[ns("selected")]] <- shiny::reactive({
            if (is.null(env$input[[ns("brush")]]))
                return( integer(0) )
        
            brush <- env$input[[ns("brush")]]
            grob <- isolate( e("grob") )
            selection <- isolate( e("selection") )
            
            n <- length(selection)
            from <- max(1, min(n, floor(brush$ymin*n+1.5)))
            to <- max(1, min(n, floor(brush$ymax*n+0.5)))
            
            if (to < from) 
                return(integer(0))
            
            original_rows <- grob$info$row_order$order[to:from]
            selection[original_rows]
        })
        
        #env$output[[ns("selected_text")]] <- shiny::renderUI({
        #    req(length(e("selected")) > 0)
        #
        #    shiny::div(
        #        shiny::h3("Selected"),
        #        shiny::pre(
        #            paste(
        #                feature_labels(env)[ e("selected") ], 
        #                collapse="\n")))
        #})
        
        # Gadget support
        env[[ns("output")]] <- shiny::reactive({
            list(
                rows_shown = e("selection"),
                rows_selected = e("selected")
            )
        })
    }

    composable_shiny_app(ui, server, output=ns("output"), title="Heatmap")
}




