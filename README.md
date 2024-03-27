# Functions to calculate the Planetary boundary status based on LPJmL
    simulation data

R package **boundaries**, version **0.1.1**

[![CRAN status](https://www.r-pkg.org/badges/version/boundaries)](https://cran.r-project.org/package=boundaries)  [![R build status](https://gitlab.pik-potsdam.de/tess/boundaries/workflows/check/badge.svg)](https://gitlab.pik-potsdam.de/tess/boundaries/actions) [![codecov](https://codecov.io/gh/tess/boundaries/branch/master/graph/badge.svg)](https://app.codecov.io/gh/tess/boundaries) 

## Purpose and Functionality

Provide calculation & visualizations functions for planetary boundaries.


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

Breier J, Stenzel F, Braun J (2024). _boundaries: Functions to calculate the Planetary boundary status based on LPJmL simulation data_. R package version 0.1.1.

A BibTeX entry for LaTeX users is

 ```latex
@Manual{,
  title = {boundaries: Functions to calculate the Planetary boundary status based on LPJmL
simulation data},
  author = {Jannes Breier and Fabian Stenzel and Johanna Braun},
  year = {2024},
  note = {R package version 0.1.1},
}
```
