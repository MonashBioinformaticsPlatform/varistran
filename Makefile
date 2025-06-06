

document :
	Rscript -e "devtools::document()"

check :
	Rscript -e "devtools::check()"

site : document
	Rscript -e "pkgdown::build_site()"

install :
	Rscript -e "devtools::install()"

paper.pdf : paper.md paper.bib
	pandoc --filter pandoc-citeproc -o paper.pdf paper.md


test : test_advice test_biplot test_heatmap test_heatmap_big test_pathological test_vst_plot test_vst_stability_plots

test_advice :
	Rscript test/advice.R

test_biplot :
	Rscript test/biplot.R

test_heatmap :
	Rscript test/heatmap.R

test_heatmap_big :
	Rscript test/heatmap_big.R

test_pathological :
	Rscript test/pathological.R

test_vst_plot :
	Rscript test/vst_plot.R

test_vst_stability_plots :
	Rscript test/vst_stability_plots.R


# This launches a web-server, so not included in standard test run
test_shiny:
	Rscript test/shiny.R










