#
# Shiny wrappers around plotting functions
#

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
        prefix=paste0(prefix,"plot")
    )

    ui <- shiny::tags$div(
        titlePanel("Stability plot"),
        numericInput(p("bins"), "Bins", value=20,min=1,max=10000),
        plot$component_ui
    )

    server <- function(env) {
        plot$component_server(env)
    }

    composable_shiny_app(ui, server)
}


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
        prefix=p("plot")
    )

    ui <- shiny::tags$div(
        titlePanel("Biplot"),
        numericInput(p("n_features"), "Number of labelled features", 20, min=0, step=1),
        numericInput(p("balance"), "Feature/sample relative scaling", 0.25, min=0,step=0.05),
        plot$component_ui
    )

    server <- function(env) {
        plot$component_server(env)
    }

    composable_shiny_app(ui, server)
}
