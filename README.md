# BDE_R_wofklow
 The R workflow to estimate global, continental, and country-level bee species richness. A full implementation of this code is available in [BeeBDC version ≥1.3.0](https://jbdorey.github.io/BeeBDC/index.html). 

 # System requirements
 Runs on Mac, Linux, or Windows and requires R (≥ 2.10) and the dependencies listed [here](https://cloud.r-project.org/web/packages/BeeBDC/index.html). See recent GitHub tests for the current version [here](https://github.com/jbdorey/BeeBDC/actions). Install time depends on the currently-installed packages but should be <10 minutes with a reasonable internet connection.

 # Instructions for use and demo
 The current version of this workflow is described within the relevant .R files. However, instructions for the generlaised-implementation are available as a [vignette here](https://jbdorey.github.io/BeeBDC/articles/speciesRichness_example.html) that also shows expected outcomes using a test dataset. The test dataset in the vignette should run in <<10 minutes. A full run of the workflow presented here, across a single computer with many cores will likely take weeks (especially at the global or continental scale with 100 replicates).
