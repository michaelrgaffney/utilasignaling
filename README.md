This repository contains all code necessary to reproduce:

[*The Natural History of Child Signals of Need in Utila, Honduras: An Exploratory Study*](https://michaelrgaffney.github.io/utilasignaling/). Michael R. Gaffney et al.

This paper is not yet published.

Instructions:

1.  Clone this repository
2.  Open the project in RStudio, Positron or `cd` into the directory and launch `R`. This will automatically bootstrap [`renv`](https://rstudio.github.io/renv/index.html).
3.  After the bootstrapping process, enter the following in the console: `renv::restore()`. This should install all the necessary packages in an isolated project-specific library.
4.  knit the `Paper.qmd` file using the RStudio GUI or Positron. This will generate the preprint file [`Paper.html`](https://michaelrgaffney.github.io/utilasignaling/), which will display in the RStudio/Positron Viewer or can be viewed in any web browser. (Note: if not using RStudio or Positron, you will need a recent version of [pandoc](https://pandoc.org) and [quarto](https://quarto.org) installed.)

Note 1: Analyses used "R version 4.4.1". You might need to install this version of R to reproduce them. [`rig`](https://github.com/r-lib/rig) is one tool for managing `R` versions.

Note 2: Some results will differ slightly from the official version because identifying information, such as child ages, has been aggregated into coarser categories, or randomized, to help preserve anonymity.
