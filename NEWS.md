# varistran 1.1.0

* Add samesum log transformation method.

# varistran 1.0.5

* plot_heatmap now has parameters baseline_to and scale_to giving control over the ranges of the scales, and parameters show_baseline and show_tree giving control over what parts of the plot are shown.

# varistran 1.0.4

* Increase robustness of heatmap to missing values.
* Baseline plot is not shown in heatmap if all zero.

# varistran 1.0.3

* No code changes. Updated README with references and supporting/contributing section.

# varistran 1.0.2

* Added n parameter to plot_heatmap, to show only the top n rows by span of expression levels.
* Remove dependency on ggdenro, which is not available in R 3.4.1.
* Added packages needed for testing to suggested dependencies: biomaRt, DESeq2, NBPSeq

