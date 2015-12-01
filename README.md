# EPVDemo

Contributors:
Dan Cervone
Alex D'Amour
Luke Bornn
Kirk Goldsberry

Demo of NBA Expected Possession Value model presented in the paper ["A Multiresolution Stochastic Process Model for Predicting NBA Possession Outcomes"](http://arxiv.org/abs/1408.0777).

The source .tex for the tutorial file EPV_demo.pdf can be built from EPV_demo.Rnw using RStudio, or the command

> Rscript -e "library(knitr); knit('./EPV_demo.Rnw')"

This operation requires the `knitr` `R` package, and takes about 1 minute on most systems. Users should be running the latest version of R with the following packages (and all of their dependencies) installed: knitr, INLA, xtable, data.table, matrixStats, RColorBrewer, fields, animation. By default, the missing packages will try to be automatically installed when compiling EPV_demo.Rnw

We would like to acknowledge STATS, LLC for allowing the public release of this data sample, and hope that this tutorial helps strengthen and broaden the impact of our paper.