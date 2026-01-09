**Materials for the 'Nonstationary spawner-recruitment dynamics for Pacific salmon workshop', February 13, 2026 - 8:30am to 3:30pm - Hyatt Regency Vancouver (Room TBA), Vancouver, BC, Canada.**

Lead: Dan Greenberg (Fisheries and Oceans Canada - Pacific)

This workshop will cover material related to assessing nonstationarity in salmon population dynamics including several exercises for participants to go through.

The basis of this workshop will rely on a package of functions to fit and visualize nonstationarity in spawner-recruit dynamics:
[samEst]([https://example.com](https://github.com/Pacific-salmon-assess/samEst))

To install samEst, we recommend using the 'remotes' package as so:

```{r} 
install.packages("remotes") 
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')
```

Note, that due to the size of the TMB models, that installation may take some time (upwards of 15-20 mins) to complete. If there are issues with installation prior to the workshop, feel free to contact Dan (dan.greenberg[at]dfo-mpo.gc.ca). Note that the package and dependencies were built on R v4.3.2 - but was tested on R v.4.5.2 for compatibility - older versions may cause issues. It is also recommended to use Rstudio, for visualizing the markdown files in HTML.

To run Bayesian models with [*Stan*](https://mc-stan.org/) (via [*rstan*](https://mc-stan.org/rstan/)), you will need to download [Rtools](https://cran.r-project.org/bin/windows/Rtools/) - of a version appropriate to your current R version - follow the instructions for installation in the link.

There will be two main exercises in the workshop, interspersed with discussion and presentations. For each exercise we include an Rmarkdown file, where participants can edit the code chunks, and alternatively a pure R script, depending on user preferences. Each folder contains the scripts for each exercise.

Exercise 1: Fitting and assessing nonstationary spawner-recruit models
  1. stable dynamics
  2. stable, but autocorrelated, dynamics
  3. non-stationary dynamics
  4. regime shifts
  5. empirical examples

Exercise 2: Time-varying reference points and management strategy evaluation
  1. time-varying Smsy and Umsy
  2. harvest control rules with 'static' vs. time-varying benchmarks
  3. assessing management outcomes and trade-offs

You can clone this repository to gain access to the materials (data, scripts, etc.) or download the folder (and unzip) via the green 'code' button above.
