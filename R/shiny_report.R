
#
# All of the plots.
#

#' Shiny component to describe the output of a call to vst.
#'
#' @param y Transformed matrix.
#'
#' @param counts Counts matrix to transform.
#'
#' @param sample_labels Sample labels.
#'
#' @param prefix Prefix for variables in environment for composable shiny app.
#'
#' @export
shiny_vst <- function(y=NULL, counts=NULL, sample_labels=NULL, prefix="") {
    ns <- shiny::NS(prefix)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)
    sample_labels <- ensure_reactable(sample_labels)

    ui <- function(request) shiny::div(
        shiny::h3("Transformation"),
        shiny::uiOutput(ns("report"))
    )

    server <- function(env) {
        env[[ns("y")]] <- shiny::reactive({
            y_val <- y(env)
            counts_val <- counts(env)
            if (is.null(y_val))
                y_val <- vst(counts_val)
            y_val
        })

        env$output[[ns("report")]] <- shiny::renderUI({
            y <- env[[ns("y")]]()

            sample_labels_val <- sample_labels(env)
            if (is.null(sample_labels_val))
                sample_labels_val <- colnames(y)
            if (is.null(sample_labels_val))
                sample_labels_val <- as.character(seq_len(ncol(y)))
            
            if (!is.null(attr(y,"method")) && attr(y,"method") == "samesum") {
                advice <- vst_advice(y)
                colnames(advice) <- c("Count","Transformed count","2-fold step")
                
                libs <- attr(y,"samples")
                colnames(libs) <- c("Sample", "Scale", "True library size", "Adjusted library size")
                
                advice_html <- list(                    
                    shiny::p("Samesum tranformation."),
                    shiny::h4("Library sizes and scaling"),
                    shiny::renderTable(libs, include.rownames=F, digits=2),
                    shiny::p(paste0("Library size adjustment method: samesum")),
                    
                    shiny::h4("Transformation"),
                    shiny::renderTable(advice, include.rownames=F, digits=2),
                    shiny::p(                        
                        "The column \"2-fold step\" shows the difference in transformed ",
                        "count from the previous row. If a simple log2 transformation were ",
                        " used this would be uniformly 1, however the addition of a pseudocount ",
                        "makes this smaller for counts close to zero."
                        ),
                    shiny::p(
                        "Note that values shown are for a typical sample."
                        ),
                    
                    parenthetically(
                        "Samesum transformation is performed ",
                        "by varistran::samesum_log2_norm or varistran::samesum_log2_cpm. ",
                        "Transformed values table produced ",
                        "by varistran::vst_advice."
                        )
                    )
            } else if (!is.null(attr(y,"method"))) {
                dispersion <- attr(y,"dispersion")
                description <- vst_methods[[attr(y,"method")]]$description
                units <- vst_methods[[attr(y,"method")]]$units
                cpm <- attr(y,"cpm")
                if (cpm)
                    units <- paste0("Units for transformed counts are ",units,
                                    " Counts Per Million.")
                else
                    units <- paste0("Units for transformed counts are ",units," read count.")
                
                libs <- data.frame(
                    Sample = sample_labels_val,
                    "True library size" = attr(y,"true.lib.size"),
                    "Adjusted library size" = attr(y,"lib.size"),
                    check.names=FALSE
                    )
            
                advice <- vst_advice(y)
                colnames(advice) <- c("Count","Transformed count","2-fold step")
                
                advice_html <- list(                    
                    shiny::p(description),
                    shiny::p(units),
                    shiny::p(sprintf("Estimated dispersion is %.4f.", dispersion)),
                    
                    shiny::h4("Library sizes"),
                    shiny::renderTable(libs, include.rownames=F, digits=0),
                    shiny::p(paste0("Library size adjustment method: ",
                        attr(y,"lib.size.method"))),
                    
                    shiny::h4("Transformation"),
                    shiny::renderTable(advice, include.rownames=F, digits=2),
                    shiny::p(                        
                        "The column \"2-fold step\" shows the difference in transformed ",
                        "count from the previous row. If a simple log2 transformation were ",
                        " used this would be uniformly 1, however a variance stabilizing ",
                        "transformation makes this smaller for counts close to zero."
                        ),
                    shiny::p(
                        "Note that values shown are for a sample with average adjusted library size sample."
                        ),
                    
                    parenthetically(
                        "Variance stabilizing transformation is performed ",
                        "by varistran::vst. ",
                        "Transformed values table produced ",
                        "by varistran::vst_advice."
                        )
                    )
            } else {
                advice_html <- list(
                    shiny::p("Any transformation not by Varistran, details unknown.")
                    )
            }

            do.call(shiny::div, advice_html)
        })
    }

    composable_shiny_app(ui, server)
}

#' Shiny component to filter rows and columns a matrix.
#'
#' @param y Transformed matrix.
#'
#' @param counts Original count matrix if available.
#'
#' @param sample_labels Sample labels.
#'
#' @param feature_labels Feature labels.
#'
#' @param prefix Prefix for variables in environment for composable shiny app.
#'
#' @export
shiny_filter <- function(y, counts=NULL, sample_labels=NULL, feature_labels=NULL, prefix="") {
    ns <- shiny::NS(prefix)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)
    sample_labels <- ensure_reactable(sample_labels)
    feature_labels <- ensure_reactable(feature_labels)

    ui <- function(request) shiny::uiOutput(ns("ui"))

    server <- function(env) {
        env[[ns("vars")]] <- shiny::reactive(shiny::withProgress(message="Loading", {
            y_val <- y(env)
            counts_val <- counts(env)
            sample_labels_val <- sample_labels(env)
            choices <- seq_len(ncol(y_val))
            if (!is.null(sample_labels_val))
                names(choices) <- sample_labels_val
            else if (!is.null(colnames(y_val)))
                names(choices) <- colnames(y_val)
            
            default_min_count <- 5
            default_min_expression <- min(y_val, 0.0)
            default_samples <- choices
            
            if (is.null(counts_val)) {
                counts_input <- ""
                what <- "expression level"
            } else {
                counts_input <- shiny::numericInput(ns("min_count"), "Minimum mean count", default_min_count)
                what <- "transformed count"
            }
            
            ui <- shiny::tags$div(
                shiny::tags$h3("Select samples"),
                shiny::selectInput(ns("samples"), "Select samples", 
                    selected=default_samples, choices=choices, multiple=TRUE),
                shiny::tags$h3("Filter features"),
                counts_input,
                shiny::numericInput(ns("min_expression"), 
                    paste0("Minimum mean ",what), default_min_expression),
                shiny::textOutput(ns("report"))
            )
            
            list(
                ui=ui,
                y=y_val,
                counts=counts_val,
                default_samples=default_samples,
                default_min_count=default_min_count,
                default_min_expression=default_min_expression
            )
        }))
        
        env$output[[ns("ui")]] <- shiny::renderUI({ env[[ns("vars")]]()$ui })

        env[[ns("filtered")]] <- shiny::reactive({
            vars <- env[[ns("vars")]]()
            y_val <- vars$y
            counts_val <- vars$counts
            
            # Might not exist if tab hasn't been viewed :-(
            samples_val <- env$input[[ns("samples")]]
            if (is.null(samples_val))
                samples_val <- vars$default_samples
            
            min_expression_val <- env$input[[ns("min_expression")]]
            if (is.null(min_expression_val))
                min_expression_val <- vars$default_min_expression
                
            min_count_val <- env$input[[ns("min_count")]]
            if (is.null(min_count_val))
                min_count_val <- vars$default_min_count
            
            sample_select <- as.numeric(samples_val)
            feature_select <- which(
                rowMeans(y_val[,sample_select,drop=FALSE]) >= min_expression_val
            )
                        
            if (!is.null(counts_val)) {
                feature_select <- intersect(feature_select,which(
                    rowMeans(counts_val[,sample_select,drop=FALSE]) >= min_count_val
                ))
            }
            
            list(
                sample_select=sample_select,
                feature_select=feature_select,
                y=y_val[feature_select,sample_select,drop=FALSE],
                counts=counts_val[feature_select,sample_select,drop=FALSE],
                sample_labels=sample_labels(env)[sample_select],
                feature_labels=feature_labels(env)[feature_select]
            )
        })

        env$output[[ns("report")]] <- shiny::renderText({
            filtered <- env[[ns("filtered")]]()
            sprintf("%d of %d features will be used.",
                length(filtered$feature_select),
                nrow(y(env))
            )
        })
    }

    composable_shiny_app(ui, server)
}


#' Shiny report
#'
#' Produce an interactive Shiny report showing diagnostic plots of transformed counts.
#'
#' @param y A matrix of exprssion levels, such as a transformed counts matrix.
#'
#' @param counts Original counts.
#'
#' @param sample_labels Optional. Sample names.
#'
#' @param feature_labels Optional. Feature names.
#'
#' @param prefix Optional, to facilitate use as a component of a larger Shiny app. Inputs and outputs are given this prefix.
#'
#' @return A shiny.appobj.
#'
#' Either y or counts or both must be given.
#'
#' Used interactively, the shiny report runs immediately. Otherwise it can be launched by print()-ing it. A call to this function can also be the last line in an app.R file in a Shiny app directory.
#'
#' @author Paul Harrison
#'
#' @examples
#'
#' # Generate some random data.
#' means <- runif(100,min=0,max=1000)
#' counts <- matrix(rnbinom(1000, size=1/0.01, mu=rep(means,10)), ncol=10)
#'
#' y <- varistran::vst(counts)
#' if (interactive())
#'     varistran::shiny_report(y, counts)
#'
#' @export
shiny_report <- function(y=NULL, counts=NULL, sample_labels=NULL, feature_labels=NULL, prefix="") {
    ns <- shiny::NS(prefix)

    y <- ensure_reactable(y)
    counts <- ensure_reactable(counts)
    sample_labels <- ensure_reactable(sample_labels)
    feature_labels <- ensure_reactable(feature_labels)

    transform <- shiny_vst(y, counts, sample_labels, prefix=ns("transform"))
    ty <- function(env) env[[ns("transform-y")]]()

    filter <- shiny_filter(ty, counts, sample_labels, feature_labels, prefix=ns("filter"))
    fy <- function(env) env[[ns("filter-filtered")]]()$y
    fcounts <- function(env) env[[ns("filter-filtered")]]()$counts
    fsample_labels <- function(env) env[[ns("filter-filtered")]]()$sample_labels
    ffeature_labels <- function(env) env[[ns("filter-filtered")]]()$feature_labels

    stability <- shiny_stability(fy, fcounts, prefix=ns("stability"))
    mds_plot <- shiny_mds_plot(fy, sample_labels=fsample_labels, prefix=ns("mds"))
    biplot <- shiny_biplot(fy, sample_labels=fsample_labels, feature_labels=ffeature_labels,  prefix=ns("biplot"))
    heatmap <- shiny_heatmap(fy, sample_labels=fsample_labels, feature_labels=ffeature_labels,  prefix=ns("heatmap"))

    panels <- list(
            function(request) shiny::tabPanel("Transform", transform$component_ui(request)),
            function(request) shiny::tabPanel("Select and filter", filter$component_ui(request)),
            function(request) shiny::tabPanel("Stability", stability$component_ui(request)),
            function(request) shiny::tabPanel("MDS plot", mds_plot$component_ui(request)),
            function(request) shiny::tabPanel("Biplot", biplot$component_ui(request)),
            function(request) shiny::tabPanel("Heatmap", heatmap$component_ui(request))
    )

    ui <- function(request) shiny::div(
        shiny::div(
            style="font-size: 150%; color: #bbbbbb; text-align: right; letter-spacing: 0.25em;",
            "Varistran"
        ),
        do.call(shiny::navlistPanel, 
            c(list(widths=c(2,10),well=FALSE), 
              lapply(panels, call_ui, request)
            )) 
    )

    server <- function(env) {
        e <- function(name) env[[ns(name)]]()
    
        transform$component_server(env)
        filter$component_server(env)
        stability$component_server(env)
        mds_plot$component_server(env)
        biplot$component_server(env)
        heatmap$component_server(env)
        
        
        env[[ns("columns")]] <- shiny::reactive({
            e("filter-filtered")$sample_select
        })
        
        env[[ns("rows")]] <- shiny::reactive({
            e("filter-filtered")$feature_select
        })
        
        env[[ns("rows-selected")]] <- shiny::reactive({
            e("rows")[ e("heatmap-selected") ]
        })
    }

    app <- composable_shiny_app(ui, server)
    app$component_panels <- panels
    
    app
}



