# Get the simulation time based on the information in the configuration file.
get_sim_time <- function(config) {
  if (is.character(config)) {
    config <- lpjmlkit::read_config(config)
  }
  return(as.character(config$outputyear:config$lastyear))
}


# Get the simulation cells based on the information in the configuration file.
get_sim_cells <- function(config, all_as_null = TRUE) {
  if (is.character(config)) {
    config <- lpjmlkit::read_config(config)
  }
  if (config$startgrid == "all" && config$endgrid == "all") {
    config$startgrid <- 0
    config$endgrid <- 67419
  }
  if (all_as_null && config$startgrid == 0 && config$endgrid == 67419) {
    return(NULL)
  }
  return(as.character(config$startgrid:config$endgrid))
}


# Estimate the spatial subset of the scenario and reference run if simulations
# are not run globally.
get_spatial_subset <- function(config_one, config_two) {

  spatial_subset_scenario <- get_sim_cells(
    config = config_one,
    all_as_null = TRUE
  )
  spatial_subset_reference <- get_sim_cells(
    config = config_two,
    all_as_null = TRUE
  )
  if (is.null(spatial_subset_scenario) && is.null(spatial_subset_reference)) {
    spatial_subset <- NULL
  } else if (is.null(spatial_subset_scenario) ||
               is.null(spatial_subset_reference)) {
    spatial_subset <- c(spatial_subset_scenario, spatial_subset_reference)
  } else {
    spatial_subset <- intersect(
      spatial_subset_scenario,
      spatial_subset_reference
    )
  }
  spatial_subset
}


# convenience function to perform spatial_subsets
conditional_subset <- function(x, spatial_subset) {
  if (is.null(spatial_subset)) {
    return(x)
  }
  return(x$subset(cell = spatial_subset))
}


# Get the simulation outputs based on the information in the configuration file.
get_sim_outputs <- function(config) {
  if (is.character(config)) {
    config <- lpjmlkit::read_config(config)
  }

  return(
    sapply( # nolint:undesirable_function_linter
      config$output,
      function(x, with_meta) {
        y <- x$file$name
        # Workaround to use the package data
        if (startsWith(y, "extdata")) {
          y <- system.file(y, package = "boundaries")
        }
        has_suffix <- grepl("\\.(nc|bin|clm)", y)
        if (with_meta && !has_suffix) {
          y <- paste0(y, config$default_suffix, ".json")
        } else if (with_meta && has_suffix) {
          y <- paste0(y, ".json")
        }
        names(y) <- x$id
        y
      },
      with_meta = config$output_metafile
    )
  )
}


# Get the simulation inputs based on the information in the configuration file.
get_sim_inputs <- function(config) {
  if (is.character(config)) {
    config <- lpjmlkit::read_config(config)
  }
  z <- list()
  for (i in seq_along(config$input)) {
    name <- config$input[[i]]$name
    if (startsWith(name, "/p/")) {
      z[[names(config$input)[i]]] <- name
      # Workaround to use the package data
    } else if (startsWith(name, "extdata")) {
      z[[names(config$input)[i]]] <- system.file(name, package = "boundaries")
    } else {
      z[[names(config$input)[i]]] <- file.path(config$inpath, name)
    }
  }
  z
}

# Get both, the simulation inputs and outputs based on the information in the
# configuration file.
get_sim_data <- function(config) {
  if (is.character(config)) {
    config <- lpjmlkit::read_config(config)
  }
  inputs <- get_sim_inputs(config)
  outputs <- get_sim_outputs(config)

  c(outputs, inputs)
}
