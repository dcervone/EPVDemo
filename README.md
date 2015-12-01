# EPVDemo

Contributors:
- Dan Cervone
- Alex D'Amour
- Luke Bornn
- Kirk Goldsberry

This repository contains data and code offering a demo of the NBA Expected Possession Value model presented in the paper ["A Multiresolution Stochastic Process Model for Predicting NBA Possession Outcomes"](http://arxiv.org/abs/1408.0777).

The main document that introduces and illustrates the code/data is EPV_demo.pdf. The source .tex for the tutorial file EPV_demo.pdf can be built from EPV_demo.Rnw using RStudio, or the command

> Rscript -e "library(knitr); knit('./EPV_demo.Rnw')"

This operation requires the `knitr` `R` package, and takes about 1 minute on most systems. Users should be running the latest version of `R` with the following packages (and all of their dependencies) installed: `sp`, `knitr`, `INLA`, `xtable`, `data.table`, `matrixStats`, `RColorBrewer`, `fields`, `animation`. By default, the missing packages will try to be automatically installed when compiling EPV_demo.Rnw

We would like to acknowledge STATS, LLC for consenting the inclusion of a full-game data sample.