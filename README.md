# <a href=''><img src='inst/img/logo.png' align='right' alt='logo' height=139 /></a> Planetary Boundary Status based on LPJmL simulations

R package **boundaries**, version **1.2.0**

[![CRAN status](https://www.r-pkg.org/badges/version/boundaries)](https://cran.r-project.org/package=boundaries) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11550559.svg)](https://doi.org/10.5281/zenodo.11550559) [![R build status](https://github.com/PIK-tess/boundaries/workflows/check/badge.svg)](https://github.com/PIK-tess/boundaries/actions) [![codecov](https://codecov.io/gh/PIK-tess/boundaries/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PIK-tess/boundaries) [![r-universe](https://pik-piam.r-universe.dev/badges/boundaries)](https://pik-piam.r-universe.dev/builds)

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

## Tutorial

The package comes with a vignette describing the basic functionality of the package and how to use it. You can load it with the following command (the package needs to be installed):

```r
vignette("calc_and_plot_boundaries") # Calculate and plot boundaries
```

## Questions / Problems

In case of questions / problems please contact Johanna Braun <braun@pik-potsdam.de>.

## Citation

To cite package **boundaries** in publications use:

Braun J, Breier J, Stenzel F, Vanelli C (2024). _boundaries: Planetary Boundary Status based on LPJmL simulations_. doi:10.5281/zenodo.11550559 <https://doi.org/10.5281/zenodo.11550559>, R package version 1.2.0, <https://github.com/PIK-tess/boundaries>.

A BibTeX entry for LaTeX users is

 ```latex
@Manual{,
  title = {boundaries: Planetary Boundary Status based on LPJmL simulations},
  author = {Johanna Braun and Jannes Breier and Fabian Stenzel and Caterina Vanelli},
  year = {2024},
  note = {R package version 1.2.0},
  url = {https://github.com/PIK-tess/boundaries},
  doi = {10.5281/zenodo.11550559},
}
```
