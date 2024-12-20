**Calculate and plot boundaries** ✎📊 describes the calculation and
plotting of planetary boundary statuses at different spatial scales
based on LPJmL. It aims to simplify the post-processing of LPJmL output
data and to make the calculations and visualizations transparent and
flexible.

## Overview

The calculations rely on the output of LPJmL simulations. To configure
and submit LPJmL simulations, the LPJmL Runner from the lpjmlkit package
can be used (<https://github.com/PIK-LPJmL/lpjmlkit>). In doing so, it
is important to ensure that all necessary output variables are written
out (at the needed spatial resolution). To do so, `boundaries` provides
a function to list all needed LPJmL input and output variables (see
`?list_outputs`). These have to be included in the LPJmL configuration
file (`lpjml_config.cjson`). Example output files can be downloaded from
<https://doi.org/10.5281/zenodo.12171096>.

For the post-processing, two working steps are required: Calculate the
boundaries (1) and plot them (2). In a seperate step, a global
validation table can be created to compare the simulated key variables
with global literature estimates (3).

#### **1. ✏ Calculate boundaries** 

  
Use `calc_status()` to calculate the boundary status. Define which
boundaries should be calculated, at which `spatial_scale` and for which
`time_span`, as well as additional arguments such as the calculation
`approach` to use (see `?calc_status` *for more information*).  
Each boundary calculation
`boundary = c("bluewater", "greenwater", "lsc", "nitrogen", "biosphere")`
has its own underlying calculation function, for which additional
arguments can be passed to `calc_status()`, e.g. the minimum tree cover
`tree_cover_thresholds` to define forest biomes in `lsc_status()`.  
See e.g. `?lsc_status` (boundary + “\_status”) for more information.  
The function returns a list with objects of class control_variable,
referring to the status of the contol variable(s) with attributes, such
as the boundary and high risk values, the unit, and the spatial
resolution.

``` r
global_status <- calc_status(
  boundary = c("bluewater", "greenwater", "lsc", "nitrogen", "biosphere"),
  config_scenario = config_path_scenario,
  config_reference = config_path_reference,
  spatial_scale = "global",
  path_baseline = output_path_baseline,
  time_span_scenario = as.character(1901:2017),
  time_span_reference = as.character(1500:1699),
  time_span_baseline = as.character(1901:2017),
  time_series_avg = 1, # computation of timeseries, no moving average
  approach = list("bluewater" = "porkka2024",
                  "greenwater" = "porkka2024",
                  "nitrogen" = "schulte_uebbing2022")
)
```

#### **2. 📊 Plot boundary status(es)** 

  
The status of the control variable(s) can be visualized using the
`plot_status()` function, by transferring the output of the
`calc_status()` function. Depending on the spatial scale, the function
applies different plotting functions (maps for `grid` and `regional`
scale; timeseries plots for `global` scale). See `?plot_status()` *for
more information*, as well as the individual plotting functions called
by `plot_status()` for the different spatial scales (`status_map()` for
`grid` and `regional` scale; `status_stylized()` (radial plot) or
`status_global` (Cartesian coordinate system) for `global` scale). For
some plots, the legend is not automatically plotted, but can be plotted
separately using the `status_legend()` function.

``` r
plot_status(global_status, filename = "./global_status.png")
status_legend(filename = "./legend.png")
```

#### **3. 📋 Create global validation table** 

  
The function `validate_simulation()` compares relevant key variables
that impact the control variable status, such as (i) inputs to the LPJmL
model (e.g. nitrogen fertilization) and (ii) processed LPJmL outputs
(e.g. irrigation water consumption), as well as (iii) the control
variable status values themselves with global literature estimates. The
parameter settings should match those chosen for application of the
`calc_status()` function at the global level, except for
`time_span_scenario`, which should match the timespans/years of the
literature estimates (see `./extdata/global_validation_data.csv` for an
overview of included literature data). The computed results will be
averaged over this timespan. See `?validate_simulation` *for more
information*.

``` r
validation_table <- validate_simulation(
  config_scenario = config_path_scenario,
  config_reference = config_path_reference,
  path_baseline = output_path_baseline,
  time_span_scenario = as.character(2005:2014),
  time_span_reference = as.character(1500:1699),
  filename = "./validation_table.csv",
  approach = list("bluewater" = "porkka2024",
                  "greenwater" = "porkka2024",
                  "nitrogen" = "schulte_uebbing2022")
)
```

#### **Miscellaneous** 

  
More helpful functions that come with the package are:

- `classify_biomes()` to determine spatially-explicit biome extents
  based on the foliar projected cover of different plant functional
  types simulated in LPJmL, the annual mean temperature and elevation
  (inputs to LPJmL). This function is used, amongst others, to calculate
  the land system change boundary status but can also be used
  independently.

- `plot_biomes()` to visualize the spatially-explicit biome extents
  calculated with `classify_biomes`.

- `calc_efrs()` to calculate environmental flow requirements based on
  different defined approaches based on LPJmL discharge output data.
  This is used for the blue water boundary status at the gridded level
  but can be relevant also for other applications.

- `as_risk_level()` to convert the status of a control variable to a
  risk level based on the boundary and high risk values. If the risk
  level is not to be plotted but needed for further analysis, this
  function can be directly applied to the output of `calc_status`.

For an overview on all implemented approaches, the respective control
variables, units, the default boundary & high risk values, as well as
all needed LPJmL outputs, see `./extdata/metric_files.yml`.

   

## Usage

``` r
library(boundaries) # nolint:undesirable_function_linter

# set paths to lpjml configuration files and outputs
path_configurations <- "./lpjml/configurations/"
config_path_scenario <- paste0(path_configurations, "config_scenario.json")
config_path_reference <- paste0(path_configurations, "config_reference.json")
path_outputs <- "./lpjml/output/"
path_outputs_pnv <- paste0(path_outputs, "pnv/")

plotpath <- paste0("./R/plots/")
```

   

### Global status timeseries 

#### **Example** *Global status calculation with moving average*

 

``` r
### calc status at the global level

# define analysis and reference period
time_span_scenario <- as.character(1901:2017)
time_span_reference <- as.character(1500:1699)

# define the number of years to average over, for a moving average; if set to
# 1, no moving average is calculated
nyear_window <- 5

global_status <- calc_status(
  # define boundaries to calculate - it can also be only one boundary
  boundary = c("bluewater", "greenwater", "lsc", "nitrogen", "biosphere"),
  config_scenario = config_path_scenario,
  config_reference = config_path_reference,
  spatial_scale = "global",
  time_span_scenario = time_span_scenario,
  time_span_reference = time_span_reference,
  # moving average over nyear_window years
  time_series_avg = nyear_window,
  approach = list("bluewater" = "porkka2024",
                  "greenwater" = "porkka2024",
                  "nitrogen" = "schulte_uebbing2022"),
  # boundary specific parameters, see individual boundary functions for details
  time_span_baseline = time_span_scenario, # for biosphere integrity
  path_baseline = paste0(path_outputs, "/pnv/"), # for biosphere integrity
  savanna_proxy = list(vegc = 7500) # for forest biome definition in the land system change boundary
)

### plot timeseries

# There are three options for plotting the global status timeseries:
# a) one timeseries panel for each boundary control variable, in a Cartesian
#    coordinate system
# b) all timeseries plotted in one panel, in a Cartesian coordinate system
# c) one timeseries panel for each boundary control variable, in a polar
#    coordinate system (called "stylized")

# a) one timeseries panel for each boundary control variable
plot_status(
  x = global_status,
  stylized = FALSE,
  filename = "./global_timeseries_panels.png"
)

# b) all timeseries plotted in one panel
plot_status(
  x = global_status,
  stylized = FALSE,
  all_in_one = TRUE,
  filename = "./global_timeseries_all_in_one.png"
)

# c) one timeseries panel for each boundary control variable, in a polar
#    coordinate system
plot_status(
  x = global_status,
  stylized = TRUE,
  filename = "./global_timeseries_stylized.png"
)
```

#### **Example** *Compare two calculation approaches for the same boundary*

 

``` r
# In this example, only the bluewater boundary is calculated, for two different
# approaches. The results are then jointly plotted in a Cartesian coordinate
# system.

# define analysis and reference period
time_span_scenario <- as.character(1901:2017)
time_span_reference <- as.character(1500:1699)

# Aprroach A: Following Porkka et al. (2024) (https://doi.org/10.1038/s44221-024-00208-7)
# Referring to global area with discharge deviations outside the pre-industrial
# range

bluewater_status_porkka <- calc_status(
  # define boundaries to calculate - it can also be only one boundary
  boundary = c("bluewater"),
  config_scenario = config_path_scenario,
  config_reference = config_path_reference,
  spatial_scale = "global",
  time_span_scenario = time_span_scenario,
  time_span_reference = time_span_reference,
  # no moving average, by setting the number of years to average over to 1:
  time_series_avg = 1,
  approach = list("bluewater" = "porkka2024")
)

# Aprroach B: Following Rockström et al. (2009) (https://doi.org/10.1038/461472a)
# Referring to global bluewater consumption, but adapting the defaullt boundary
# and high risk values

bluewater_status_rockstroem <- calc_status(
  # define boundaries to calculate - it can also be only one boundary
  boundary = c("bluewater"),
  config_scenario = config_path_scenario,
  config_reference = config_path_reference,
  spatial_scale = "global",
  time_span_scenario = time_span_scenario,
  time_span_reference = time_span_reference,
  # no moving average, by setting the number of years to average over to 1:
  time_series_avg = 1,
  approach = list("bluewater" = "rockstroem2009"),
  # change default planetary boundary and high risk values, following Gerten
  # et al. 2013 (https://doi.org/10.1016/j.cosust.2013.11.001),
  # all referring to km3/yr as defined in `metric_files.yml`
  thresholds = list("bluewater" = list(holocene = 0,
                                       pb = 2800,
                                       high_risk = 4000))
)

# plot both results in one timeseries panel
# merge both results into one list
# TODO test if this is working!
bluewater_status <- list("bluewater" = bluewater_status_porkka$bluewater,
                         "bluewater" = bluewater_status_rockstroem$bluewater)

plot_status(
  x = bluewater_status,
  stylized = FALSE,
  all_in_one = TRUE,
  filename = "./global_bluewater_status_porkka_vs_rockstroem.png"
)
```

### Gridded status 

#### **Example** *Gridded status calculation and plotting*

 

``` r


# define analysis and reference period
time_span_scenario <- as.character(2008:2017)
time_span_reference <- as.character(1500:1699)
# if the land use effect is to be isolated, the time span for the reference
# period can be set to the same as the scenario period:
# time_span_reference <- as.character(2008:2017) #nolint
# This way, climate change induced changes are excluded.

# calc status at the gridded level
gridded_status <- calc_status(
  boundary = c("lsc", "nitrogen", "bluewater", "biosphere"),
  config_scenario = config_path_scenario,
  config_reference = config_path_reference,
  time_span_scenario = time_span_scenario,
  time_span_reference = time_span_reference,
  spatial_scale = "grid",
  # set time_series_avg to NULL, to not calculate a timeseries, but to average
  # over the entire scenario time span
  time_series_avg = NULL,
  # boundary specific parameters, see individual boundary functions for details
  time_span_baseline = time_span_scenario, # for biosphere integrity
  path_baseline = paste0(path_outputs, "/pnv/"), # for biosphere integrity
  savanna_proxy = list(vegc = 7500) # for forest biome definition in the land system change boundary
)

# plot status at the gridded level, for details see `?status_maps`
# There are two options for plotting a map with the gridded status(es):
# a) plot the control variable status of each boundary (e.g. deforestation share
#   for the land system change boundary)
# b) plot the risk level of each boundary, based on a normalized color scale
#    and the boundary and high risk values

# a) control variable status
plot_status(
  x = gridded_status,
  filename = "./gridded_status_control_variable.png",
  grid_path = paste0(path_outputs_pnv, "grid.bin.json"),
  risk_level = FALSE
)

# b) risk level
plot_status(
  x = gridded_status,
  filename = "./gridded_status_risk_level.png",
  grid_path = paste0(path_outputs_pnv, "grid.bin.json"),
  risk_level = TRUE
)
```

     

### Regional status 

For the regional status, the calculation is performed on a aggregated
regional level that makes sense for the respective boundary. For
example, for the land system change boundary, the regional status is
calculated at the level of forest biomes, for biosphere integrity at the
level of biomes, and for bluewater at the level of river basins.

#### **Example** *Regional status calculation and plotting*

 

``` r

# define analysis and reference period
time_span_scenario <- as.character(2008:2017)
time_span_reference <- as.character(1500:1699)

# calc status at the regional level
regional_status <- calc_status(
  boundary = c("lsc", "bluewater", "greenwater", "biosphere"),
  config_scenario = config_path_scenario,
  config_reference = config_path_scenario,
  time_span_scenario = time_span_scenario,
  time_span_reference = time_span_reference,
  spatial_scale = "regional",
  approach = list("bluewater" = "porkka2024",
                  "greenwater" = "porkka2024"),
  path_baseline = paste0(path_outputs, "/pnv/"), # for biosphere integrity
  time_span_baseline = time_span_scenario, # for biosphere integrity
  savanna_proxy = list(vegc = 7500), # for forest biome definition in the land system change boundary
)

# plot status at the regional level

# As for the gridded status, there are two options for plotting a map with the
# regional status(es): directly the control variable status or the risk level.

# risk level status
plot_status(
  x = regional_status,
  filename = "./regional_status_risk_level.png",
  grid_path = paste0(path_outputs_pnv, "grid.bin.json"),
  risk_level = TRUE
)
```

   

### Notes & tips

- Some boundary status calculations may take long (particularly for
  green and bluewater based on the approach following Porkka et
  al. 2024). To speed up the calculation, the `in_parallel` argument can
  be set to `TRUE` in the `calc_status` function to parallelize some of
  the calculations.
- It is advisable to submit the calculations as a job to slurm
