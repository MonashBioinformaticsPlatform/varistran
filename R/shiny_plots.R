#
# Shiny wrappers around plotting functions
#

#' @export
shiny_stability <- function(y, x=NULL, design=NULL, bins=20, prefix="") {
    p <- function(name) paste0(prefix,name)

    y <- ensure_reactable(y)
    x <- ensure_reactable(x)
    design <- ensure_reactable(design)
    bins <- ensure_reactable(bins)

    plot <- shiny_plot(
        function(env) {
            bins <- env$input[[p("bins")]]
            stopifnot(bins >= 1)
            print(plot_stability(
                y=y(env),
                x=x(env),
                design=design(env),
                bins=bins
            ))
        },
        dlname="stability_plot",
        prefix=p("plot_")
    )

    ui <- shiny::tags$div(
        shiny::titlePanel("Stability plot"),
        shiny::numericInput(p("bins"), "Bins", value=20,min=1,max=10000),
        plot$component_ui
    )

    server <- function(env) {
        plot$component_server(env)
    }

    composable_shiny_app(ui, server)
}


#' @export
shiny_biplot <- function(x, sample_labels=NULL, feature_labels=NULL, n_features=20, balance=0.25, text_size=0.025, prefix="") {
    p <- function(name) paste0(prefix,name)

    x <- ensure_reactable(x)
    sample_labels <- ensure_reactable(sample_labels)
    feature_labels <- ensure_reactable(feature_labels)

    plot <- shiny_plot(
        function(env) {
            print(plot_biplot(
                x = x(env),
                sample_labels = sample_labels(env),
                feature_labels = feature_labels(env),
                n_features = env$input[[p("n_features")]],
                balance = env$input[[p("balance")]],
                text_size = 0.025
            ))
        },
        dlname="biplot",
        prefix=p("plot_")
    )

    ui <- shiny::tags$div(
        shiny::titlePanel("Biplot"),
        shiny::numericInput(p("n_features"), "Number of labelled features", 20, min=0, step=1),
        shiny::numericInput(p("balance"), "Feature/sample relative scaling", 0.25, min=0,step=0.05),
        plot$component_ui
    )

    server <- function(env) {
        plot$component_server(env)
    }

    composable_shiny_app(ui, server)
}

#' @export
shiny_heatmap <- function(y, sample_labels=NULL, feature_labels=NULL, units="units", prefix="") {
    p <- function(name) paste0(prefix,name)

    y <- ensure_reactable(y)
    sample_labels <- ensure_reactable(sample_labels)
    feature_labels <- ensure_reactable(feature_labels)
    units <- ensure_reactable(units)

    plot <- shiny_plot(
        callback = function(env) {
            print(env[[p("grob")]]())
        },
        width=500,
        height=700,
        dlname="heatmap",
        prefix=p("plot_")
    )

    ui <- shiny::tags$div(
        shiny::titlePanel("Heatmap"),
        "Features are selected based on variance around their mean expression level.",
        shiny::numericInput(p("n"), "Numer of features to show", 50, min=10,max=2000,step=10),
        shiny::checkboxInput(p("cluster_samples"), "Cluster samples", FALSE),
        plot$component_ui
    )

    server <- function(env) {
        env[[p("grob")]] <- reactive({
            y_val <- as.matrix(y(env))
            y_centered <- y_val - rowMeans(y_val)
            y_var <- rowMeans(y_centered*y_centered)
            selection <- rep(FALSE,nrow(y_val))
            selection[ order(-y_var)[ seq_len(env$input[[p("n")]]) ] ] <- TRUE
            plot_heatmap(
                y=y_val[selection,,drop=FALSE],
                sample_labels=sample_labels(env)[selection],
                feature_labels=feature_labels(env)[selection],
                cluster_samples=env$input[[p("cluster_samples")]],
                units=units(env)
            )
        })

        plot$component_server(env)
    }

    composable_shiny_app(ui, server)
}
