# Planetary Boundary Status with LPJmL simulations <a href=''><img src='inst/img/logo.png' align='right' height='139' /></a>

R package **boundaries**, version **0.2.1**

[![CRAN status](https://www.r-pkg.org/badges/version/boundaries)](https://cran.r-project.org/package=boundaries)  [![R build status](https://gitlab.pik-potsdam.de/tess/boundaries/workflows/check/badge.svg)](https://gitlab.pik-potsdam.de/tess/boundaries/actions) [![codecov](https://codecov.io/gh/tess/boundaries/branch/master/graph/badge.svg)](https://app.codecov.io/gh/tess/boundaries) 

## Purpose and Functionality

A systematic approach to quantify the status of the terrestrial
    planetary boundaries based on the Dynamic Global Vegetation Model (DGVM)
    Lund-Potsdam-Jena managed Land (LPJmL) hosted at the Potsdam Institute for
    Climate Impact Research (PIK). The supported planetary boundaries are
    "biosphere integrity", "land-system change", "bluewater", "greenwater" and
    "nitrogen flows".
##  Overview

### &#127758;&#127777;  Terrestrial Planetary Boundary Status
`calc_status()` calculate the status of the boundaries based on &#127793; **LPJmL Simulations** depending on ...
1. &#127899; Scenario configuration
2. &#127760; Spatial scale: global, sub-global, grid
3. &#128346; Time span and resolution
4. &#128209; Methodological approach

and returns the status of each underlying control variable.


### &#127912; Status Visualization
`plot_status()` visualize the &#127777; status of the &#127758; boundaries based on the calculated control variables either as a &#128506; map or a &#128200; time series plot.


## Installation

For installation of the most recent package version an additional repository has to be added in R:

```r
options(repos = c(CRAN = "@CRAN@", pik = "https://rse.pik-potsdam.de/r/packages"))
```
The additional repository can be made available permanently by adding the line above to a file called `.Rprofile` stored in the home folder of your system (`Sys.glob("~")` in R returns the home directory).

After that the most recent version of the package can be installed using `install.packages`:

```r 
install.packages("boundaries")
```

Package updates can be installed using `update.packages` (make sure that the additional repository has been added before running that command):

```r 
update.packages()
```

## Questions / Problems

In case of questions / problems please contact Jannes Breier <jannesbr@pik-potsdam.de>.

## Citation

To cite package **boundaries** in publications use:

Breier J, Stenzel F, Braun J (2024). _boundaries: Planetary Boundary Status with LPJmL simulations_. R package version 0.2.1.

A BibTeX entry for LaTeX users is

 ```latex
@Manual{,
  title = {boundaries: Planetary Boundary Status with LPJmL simulations},
  author = {Jannes Breier and Fabian Stenzel and Johanna Braun},
  year = {2024},
  note = {R package version 0.2.1},
}
```
