#' Plot polar boundaries plot including time series of boundaries
#'
#' Plot time series of boundaries into iconic polar boundaries plot only
#' focussing on the terrestrial boundaries (half-circle). Wedges are
#' scaled and normalized according to each boundary.
#'
#' @param x list with global output from calc_status
#'
#' @param filename character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen)
#'
#' @param add_legend logical, specify whether a legend should be plotted
#'
#' @param normalization see [`as_risk_level()`]
#'
#' @examples
#' \dontrun{
#' pb_status <- calc_status(
#'   boundary = c("lsc", "biosphere", "bluewater", "greenwater", "nitrogen"),
#'   config_scenario = "./config_lu_1500_2016.json",
#'   config_reference = "./config_pnv_1500_2016.json",
#'   time_span_scenario =  as.character(1986:2016),
#'   time_span_reference =  as.character(1986:2016),
#'   spatial_scale = "global",
#'   method = list(bluewater = "porkka2023",
#'                 nitrogen = "schulte_uebbing2022"),
#'   savanna_proxy = list(vegc = 7500),
#'   time_aggregation_args = 1,
#'   path_baseline = "./pnv_1500_2016/",
#' )
#'
#' ggplot2::ggsave(
#'   "test_boundaries_plot.png",
#'   plot_status_stylized(pb_status),
#'   width = 16,
#'   height = 9,
#'   dpi = 300
#' )
#' }
#'
#' @md
#' @export
plot_status_stylized <- function(
  x,
  filename,
  add_legend = TRUE,
  normalization = "increasing risk",
  high_risk = 3.5,
  background_alpha = 1
) {

  if (add_legend) {
    ggplot2::ggsave(
      filename,
      cowplot::ggdraw(plot_legend()) + cowplot::draw_plot(
        draw_stylized(
          x,
          normalization = normalization,
          high_risk = high_risk,
          background_alpha = background_alpha
        ),
        vjust = -0.1,
        hjust = 0.025,
        scale = 0.85
      ),
      width = 16,
      height = 9,
      dpi = 300
    )
  } else {
    ggplot2::ggsave(
      filename,
      cowplot::ggdraw() + cowplot::draw_plot(
        draw_stylized(
          x,
          normalization = normalization,
          high_risk = high_risk,
          background_alpha = background_alpha
        ),
        scale = 0.85
      ),
      width = 16,
      height = 9,
      dpi = 300
    )
  }
}


draw_stylized <- function(
  x,
  normalization = "increasing risk",
  high_risk = 3.5,
  background_alpha = 1
) {

  # Convert control variable to risk_level and normalize for plotting
  x_lvl <- as_risk_level(
    x,
    type = "continuous",
    normalize = normalization
  )

  # quick fix to not allow na, negative, or values > 3
  x_lvl <- lapply(x_lvl, function(x) {
    x[x < 0 | is.nan(x) | is.na(x)] <- 0
    x
  })

  # get max y value for plotting + 5% margin
  max_y <- x_lvl %>%
    unlist() %>%
    max() %>%
    sum(0.05 * .)

  # Create tibble for plotting
  x_table <- tibble::tibble(
    x = seq_along(x_lvl),
    y = rep(max_y, length(x_lvl)),
    start_year = min(as.integer(names(x_lvl[[names(x_lvl)[1]]]))),
    end_year = max(as.integer(names(x_lvl[[names(x_lvl)[1]]]))),
    name = names(x_lvl),
    status = x_lvl
  )

  # get year in which boundary is transgressed
  x_table$transgression_year <- sapply(
    x_table$name,
    function(x) {
      as.integer(
        names(which(x_lvl[[x]] > attributes(x_lvl[[x]])$thresholds$pb)[1])
      )
    }
  )

  # if boundary is not transgressed or near to start_year end_year set alpha to
  #   1 so label for transgression year is not drawn (else 0)
  x_table$alpha <- sapply(
    x_table$transgression_year,
    function(x) {
      ifelse(
        abs(x_table$start_year[1] - x) < 20 || abs(x - x_table$end_year[1]) < 20, # nolint
        0,
        1
      )
    }
  )

  # Create ggproto class to plot time series of boundaries into wedge
  TSBoundaryStatus <- ggplot2::ggproto( # nolint
    "TSBoundaryStatus",
    ggplot2::Stat,
    compute_group = function(data, scales, span_x = 1) {
      wedge_frame <- data.frame(
        x_map = c(
          0.05, 0.95,
          rev(seq(0.05, 0.95, by = (0.95 - 0.05) / (length(data$status[[1]]) - 1))) # nolint
        ),
        y = c(
          0, 0,
          rep(1, length(data$status[[1]]))
        )
      )
      new_x <- (data$x - span_x / 2) + (span_x * wedge_frame$x_map)
      new_y <- c(data$y * wedge_frame$y[1:2], rev(data$status[[1]]))
      new_wedge <- data.frame(x = new_x, y = new_y)
      new_wedge
    },
    required_aes = c("x", "y", "status")
  )

  # Create ggproto class to plot negative/complementary part of wedge for
  #   boundaries
  TSBoundaryStatusNegative <- ggplot2::ggproto( # nolint
    "TSBoundaryStatusNegative",
    ggplot2::Stat,
    compute_group = function(data, scales, span_x = 1) {
      wedge_frame <- data.frame(
        x_map = c(
          rev(seq(0.05, 0.95, by = (0.95 - 0.05) / (length(data$status[[1]]) - 1))), #nolint
          0.05, 0.95
        ),
        y = c(
          1.05, 1.05,
          rep(1, length(data$status[[1]]))
        )
      )
      new_x <- (data$x - span_x / 2) + (span_x * wedge_frame$x_map)
      new_y <- c(rev(data$status[[1]]), data$y * wedge_frame$y[1:2])
      new_wedge <- data.frame(x = new_x, y = new_y)
      new_wedge
    },
    required_aes = c("x", "y", "status")
  )

  # Create geom function for ggplot to draw time series of boundaries into wedge
  geom_boundaries <- function(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    rule = "evenodd",
    ...,
    na.rm = FALSE, # nolint
    show.legend = NA, # nolint
    inherit.aes = TRUE, # nolint
    negative = FALSE
  ) {
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = if (negative) TSBoundaryStatusNegative else TSBoundaryStatus, # nolint
      geom = ggplot2::GeomPolygon,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm,
        rule = rule,
        ...
      )
    )
  }

  # Create table with distribution of safe space values
  safe_space <- tibble::tibble(x = 1:5, y = rep(1, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(0.3, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  # Create table with distribution of values below safe space (white circle)
  no_space <- tibble::tibble(x = 3, y = 0.29) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(0, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x + 2.6, x = x - 2.6, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  # Create table with distribution of high risk values
  increasing_risk <- tibble::tibble(x = 1:5, y = rep(high_risk, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(1, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  if (max_y >= high_risk) {
    # Create table with distribution of very high risk values
    ultra_risk <- tibble::tibble(x = 1:5, y = rep(max_y, 5)) %>%
      dplyr::mutate(vals = purrr::map(y, ~seq(high_risk, .x, by = 0.01))) %>%
      tidyr::unnest(cols = c(vals)) %>%
      dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
      dplyr::select(-y) %>%
      dplyr::rename(y = vals)
  } else {
    # else create table with distribution of very high risk values larger than
    #   maximum values to be covered by white overlay
    ultra_risk <- tibble::tibble(x = 1:5, y = rep(high_risk + high_risk * 0.1, 5)) %>%
      dplyr::mutate(vals = purrr::map(y, ~seq(max_y, .x, by = 0.01))) %>%
      tidyr::unnest(cols = c(vals)) %>%
      dplyr::mutate(xend = x + 0.45, x = x - 0.45, yend = vals) %>%
      dplyr::select(-y) %>%
      dplyr::rename(y = vals)
  }

  # Create ggplot, first with the segments that include the colour gradient for
  #   each level of risk (safe space, uncertainty zone, high risk)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = increasing_risk,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
      size = 2
    ) +
    # Use inferno colour scale for colours beyond safe space
    ggplot2::scale_colour_viridis_c(
      option = "inferno",
      begin = 0,
      end = 1,
      direction = -1
    )

  if (max_y >= high_risk) {
    # Reset colour scale for next segment to be filled with new colour gradient
    p <- p + ggnewscale::new_scale_color() +
      ggplot2::geom_segment(
        data = ultra_risk,
        ggplot2::aes(
          x = x, xend = xend, y = y, yend = yend, color = y
        ),
        size = 2
      ) +
      # Use high inferno colours for colours beyond high risk
      ggplot2::scale_colour_viridis_c(
        option = "inferno",
        begin = 0,
        end = 0,
        direction = -1
      )
  } else {
    p <- p + ggplot2::geom_segment(
      data = ultra_risk,
      ggplot2::aes(
        x = x, xend = xend, y = y + 0.02, yend = yend + 0.02, color = y
      ),
      size = 1,
      color = "white"
    )
  }

  # Reset colour scale for next segment to be filled with new colour gradient
  p <- p + ggnewscale::new_scale_color() +
    ggplot2::geom_segment(
      data = safe_space,
      ggplot2::aes(
        x = x, xend = xend, y = y - 0.01, yend = yend - 0.01, color = y
      ),
      size = 2
    ) +
    # Use green colour scale for safe space
    ggplot2::scale_color_gradient2(
      low = green,
      mid = green,
      high = green,
      midpoint = max(safe_space$y) / 7.5
    ) +
    # Hide segment area outside of boundaries with white polygon
    #   (complementary part of wedge)
    geom_boundaries(
      data = x_table,
      ggplot2::aes(
        x = x,
        y = ifelse(y >= max_y, y, max_y),
        group = x,
        status = status
      ),
      fill = "white",
      alpha = background_alpha,
      negative = TRUE
    ) +
    # Draw outline around actual boundaries time series wedge
    geom_boundaries(
      data = x_table,
      ggplot2::aes(x = x, y = y, group = x, status = status),
      color = "black",
      size = 0.6,
      fill = "transparent"
    ) +
    # Limit x-axis to a half circle
    ggplot2::scale_x_continuous(
      limits = c(0, 10)
    ) +
    # Add polar grid
    ggplot2::coord_polar(start = pi * 1.4) +
    # Remove axis labels, grid and ticks
    ggplot2::theme_void() +
    # Remove legend
    ggplot2::theme(legend.position = "none")


  # Add lines/segments to indicate the different levels of risk
  p <- p +
    # add glow effect to the safe space line
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 1, yend = 1),
      size = 1.3,
      color = darkgreen,
      alpha = 0.8,  # Add transparency
      linetype = 1
    ) +
    # Add shading for safe space line
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 0.983, yend = 0.983),
      size = 1,
      color = "white",
      alpha = 0.3,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 0.987, yend = 0.987),
      size = 0.55,
      color = "white",
      alpha = 0.6,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 0.9885, yend = 0.9885),
      size = 0.2,
      color = "white",
      alpha = 0.9,  # Add transparency
      linetype = 1
    ) +
    # add safe space line
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 1, yend = 1),
      size = 1,
      color = darkgreen,
      linetype = 1
    ) +
    # add high risk line
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 2, yend = 2),
      size = 0.6,
      color = orange,
      alpha = 0.4,  # Add transparency
      linetype = 1
    ) +
    # add line for max y value for plotting + 5% margin
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = max_y, yend = max_y),
      size = 0.8,
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      linetype = 1
    )

  # Add tickmarks for start_year, transgression_year and end_year
  p <- p +
    # start_year tickmark
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x - 0.45, y = 0, yend = max_y + 0.12),
      size = 0.8,
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      linetype = 1
    ) +
    # transgression_year tickmark
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(
        x = x + (0.43 - ((0.43 * 2) * (end_year - transgression_year) / (end_year - start_year))), # nolint
        xend = x + (0.43 - ((0.43 * 2) * (end_year - transgression_year) / (end_year - start_year))), # nolint
        y = max_y, yend = max_y + 0.12,
        alpha = alpha
      ),
      size = 1,
      # alpha = 0.7,  # Add transparency
      linetype = 1,
      color = darkgreen
    ) +
    # end_year tickmark
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x + 0.45, xend = x + 0.45, y = 0, yend = max_y + 0.12),
      size = 0.8,
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      linetype = 1
    )

  # Add lines/segments to seperate the different boundaries
  p <- p + ggplot2::geom_segment(
    data = data.frame(x = seq(0.5, 5.5, 1)),
    ggplot2::aes(x = x, xend = x, y = 0, yend = max_y + 0.5),
    size = 1.2,
    color = "black",
    linetype = 1
  ) +
    ggplot2::geom_segment(
      data = no_space,
      ggplot2::aes(
        x = x, xend = xend, y = y, yend = yend, color = y
      ),
      size = 3.5,
      color = "white"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0.44, xend = 5.56, y = 0.33, yend = 0.33),
      size = 1.2,
      color = "black",
      linetype = 1
    )

  # Add texts on the plot
  p <- p +
    # Add start_year label
    ggplot2::geom_text(
      data = x_table,
      ggplot2::aes(
        x = x - 0.41,
        y = max_y + ifelse(max_y >= 3, 0.3, 0.25),
        label = start_year
      ),
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    # Add transgression_year label
    ggplot2::geom_text(
      data = x_table,
      ggplot2::aes(
        x = x + (0.43 - ((0.43 * 2) * (end_year - transgression_year) / (end_year - start_year))), # nolint
        y = max_y + ifelse(max_y >= 3, 0.3, 0.25),
        label = transgression_year,
        alpha = alpha
      ),
      fontface = "bold",
      color = darkgreen,
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    # Add end_year label
    ggplot2::geom_text(
      data = x_table,
      ggplot2::aes(
        x = x + 0.41,
        y = max_y + ifelse(max_y >= 3, 0.3, 0.25),
        label = end_year
      ),
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    # Add boundary label on the plot
    ggplot2::geom_label(
      data = x_table,
      ggplot2::aes(
        x = x,
        y = max_y + ifelse(max_y >= 3, 1.2, 0.8),
        label = name
      ),
      color = "black",
      size = 6  # Increase the size of the labels
    )

  # Reformat ggplot for better half circle display
  #   install.packages("remotes")
  #   remotes::install_github("yjunechoe/ggtrace")
  final <- ggtrace::with_ggtrace(
    x = p + ggplot2::theme(aspect.ratio = .52),
    method = ggplot2::Layout$render,
    trace_steps = 5L,
    trace_expr = quote({
      panels <- lapply(
        panels,
        grid::editGrob,
        vp = grid::viewport(yscale = c(0.48, 1)))
    }),
    out = "g"
  )
  return(final)
}
