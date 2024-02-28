#' Plot polar boundaries plot including time series of boundaries
#'
#' Plot time series of boundaries into iconic polar boundaries plot only
#' focussing on the terrestrial boundaries (half-circle). Wedges are
#' scaled and normalized according to each boundary.
#'
#' @param x boundaries object
#'
#' @param normalization see [`as_risk_level()`]
#'
#' @examples
#' \dontrun{
#' data_path_ma <- (
#'   "/p/projects/open/Johanna/boundaries/R/r_out/r_data/global_timeseries_avg_10_porkka.RData"
#' )
#'
#' data <- load(data_path_ma) |>
#'   get()
#'
#' ggplot2::ggsave(
#'   "test_boundaries_plot.png",
#'   plot_boundaries(data),
#'   width = 16,
#'   height = 9,
#'   dpi = 300
#' )
#' }
#'
#' @md
#' @export

plot_boundaries <- function(
  x,
  normalization = "increasing risk",
  add_legend = TRUE
) {

  x_lvl <- x
  # Convert control variable to risk_level and normalize for plotting
  for (i in seq_len(length(x_lvl))) {
    x_lvl[[i]] <- as_risk_level(
      x_lvl[[i]],
      type = "continuous",
      normalize = normalization
    )
  }

  # quick fix to not allow na, negative, or values > 3
  x_lvl <- lapply(x_lvl, function(x) {
    x[x < 0 | is.nan(x) | is.na(x)] <- 0
    x
  })

  max_y <- x_lvl |>
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

  x_table$mid_year <- sapply(
    x_table$name,
    function(x) {
      as.integer(
        names(which(x_lvl[[x]] > attributes(x_lvl[[x]])$thresholds$pb)[1])
      )
    }
  )

  x_table$alpha <- sapply(
    x_table$mid_year,
    function(x) {
      ifelse(
        abs(x_table$start_year[1] - x) < 20 || abs(x - x_table$end_year[1]) < 20,
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
    dplyr::mutate(vals = purrr::map(y, ~seq(0, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  # Create table with distribution of high risk values
  increasing_risk_1 <- tibble::tibble(x = 1:5, y = rep(2.5, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(1, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  # Create table with distribution of high risk values
  increasing_risk_2 <- tibble::tibble(x = 1:5, y = rep(3, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(2.5, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  if (max_y >= 3) {
    # Create table with distribution of high risk values
    high_risk <- tibble::tibble(x = 1:5, y = rep(3.5, 5)) %>%
      dplyr::mutate(vals = purrr::map(y, ~seq(3, .x, by = 0.01))) %>%
      tidyr::unnest(cols = c(vals)) %>%
      dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
      dplyr::select(-y) %>%
      dplyr::rename(y = vals)

    if (max_y >= 3.5) {
      # Create table with distribution of high risk values
      ultra_risk <- tibble::tibble(x = 1:5, y = rep(max_y, 5)) %>%
        dplyr::mutate(vals = purrr::map(y, ~seq(3.5, .x, by = 0.01))) %>%
        tidyr::unnest(cols = c(vals)) %>%
        dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
        dplyr::select(-y) %>%
        dplyr::rename(y = vals)
    }
  }

  # Create ggplot, first with the segments that include the colour gradient for
  #   each level of risk (safe space, uncertainty zone, high risk)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = increasing_risk_1,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
      size = 2
    ) +
    ggplot2::scale_color_gradient2(
      low = yellow, mid = yellow, high = red,
      midpoint = 1.5
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_segment(
      data = increasing_risk_2,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
      size = 2
    ) +
    ggplot2::scale_color_gradient2(
      low = red, mid = red, high = purple,
      midpoint = 2.5
    )

    if (max_y >= 3) {
      # Reset colour scale for next segment to be filled with new colour gradient
      p <- p + ggnewscale::new_scale_color() +
        ggplot2::geom_segment(
          data = high_risk,
          ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
          size = 2
        ) +
        ggplot2::scale_color_gradient2(
          low = purple, mid = darkpurple, high = darkpurple,
          midpoint = 3.5
        )
      if(max_y >= 3.5) {
        # Reset colour scale for next segment to be filled with new colour gradient
        p <- p + ggnewscale::new_scale_color() +
          ggplot2::geom_segment(
            data = ultra_risk,
            ggplot2::aes(
              x = x, xend = xend, y = y, yend = yend, color = y),
            size = 2
          ) +
          ggplot2::scale_color_gradient(
            low = darkpurple, high = darkpurple
          )
      }
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
    ggplot2::scale_color_gradient2(
      low = "white",
      mid = "white",
      high = green,
      midpoint = max(safe_space$y) / 7.5
    ) +
    # Hide segment area outside of boundaries with white polygon
    #   (complementary part of wedge)
    geom_boundaries(
      data = x_table,
      ggplot2::aes(
        x = x,
        y = ifelse(y > 3.5, y, ifelse(y > 3, 3.5, 3)),
        group = x,
        status = status
      ),
      fill = "white",
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
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 1, yend = 1),
      size = 1.3,
      color = darkgreen,
      alpha = 0.8,  # Add transparency
      linetype = 1
    ) +
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
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 1, yend = 1),
      size = 1,
      color = darkgreen,
      linetype = 1
    ) +
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = 2, yend = 2),
      size = 0.6,
      color = orange,
      alpha = 0.4,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = max_y, yend = max_y),
      size = 0.8,
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      linetype = 1
    )

  p <- p +
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x - 0.45, xend = x - 0.45, y = 0, yend = max_y + 0.12),
      size = 0.8,
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(
        x = x + (0.43 - ((0.43 * 2) * (end_year - mid_year) / (end_year - start_year))),
        xend = x + (0.43 - ((0.43 * 2) * (end_year - mid_year) / (end_year - start_year))),
        y = max_y, yend = max_y + 0.12,
        alpha = alpha
      ),
      size = 1,
      # alpha = 0.7,  # Add transparency
      linetype = 1,
      color = darkgreen
    ) +
    ggplot2::geom_segment(
      data = x_table,
      ggplot2::aes(x = x + 0.45, xend = x + 0.45, y = 0, yend = max_y + 0.12),
      size = 0.8,
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      linetype = 1
    )

  # Add lines/segments to seperate the different boundaries
  p <- p +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0.5, xend = 0.5, y = 0, yend = max_y + 0.5),
      size = 1.2,
      color = "black",
      linetype = 1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 1.5, xend = 1.5, y = 0, yend = max_y + 0.5),
      size = 1.2,
      color = "black",
      linetype = 1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 2.5, xend = 2.5, y = 0, yend = max_y + 0.5),
      size = 1.2,
      color = "black",
      linetype = 1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 3.5, xend = 3.5, y = 0, yend = max_y + 0.5),
      size = 1.2,
      color = "black",
      linetype = 1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 4.5, xend = 4.5, y = 0, yend = max_y + 0.5),
      size = 1.2,
      color = "black",
      linetype = 1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 5.5, xend = 5.5, y = 0, yend = max_y + 0.5),
      size = 1.2,
      color = "black",
      linetype = 1
    )

  # Add texts on the plot
  p <- p +
    ggplot2::geom_text(
      data = x_table,
      ggplot2::aes(
        x = x - 0.43,
        y = max_y + ifelse(max_y >= 3, 0.3, 0.25),
        label = start_year
      ),
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      data = x_table,
      ggplot2::aes(
        x = x + (0.43 - ((0.43 * 2) * (end_year - mid_year) / (end_year - start_year))),
        y = max_y + ifelse(max_y >= 3, 0.3, 0.25),
        label = mid_year,
        alpha = alpha
      ),
      fontface = "bold",
      color = darkgreen,
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      data = x_table,
      ggplot2::aes(
        x = x + 0.43,
        y = max_y + ifelse(max_y >= 3, 0.3, 0.25),
        label = end_year
      ),
      color = "#464646",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    )

  # Add texts on the plot
  p <- p +
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



# Colours adapted from PIK inhouse PB graphics
purple <- rgb(117. / 255., 24. / 255., 81. / 255., 1)
darkpurple <- rgb(44. / 255., 10. / 255., 34. / 255., 1)
lightpurple <- rgb(164. / 255., 133. / 255., 153. / 255., 1)
red <- rgb(197. / 255., 18. / 255., 22. / 255., 1)
lightred <- rgb(244. / 255., 157. / 255., 136. / 255., 1)
lightgreen <- rgb(184. / 255., 213. / 255., 131. / 255., 1)
darkgreen <- rgb(84. / 255., 156. / 255., 3. / 255., 1)
green <- rgb(134. / 255., 189. / 255., 36. / 255., 1)
yellow <- rgb(251. / 255., 206. / 255., 0. / 255., 1)
orange <- rgb(0.8763552479815456, 0.4319876970396003, 0.04398308342945019, 1.0)


plot_legend <- function() {

  # Create table with distribution of safe space values
  safe_space <- tibble::tibble(x = 7.5, y = rep(3, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(2, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x, x = x - 0.13, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  safe_space_arrow <- safe_space
  safe_space_arrow$x <- 7.55
  safe_space_arrow$xend <- 7.56
  safe_space_arrow <- dplyr::filter(safe_space_arrow, y < 2.4)

  # Create table with distribution of high risk values
  increasing_risk_1 <- tibble::tibble(x = 7.5, y = rep(4.5, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(3, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x, x = x - 0.13, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  # Create table with distribution of high risk values
  increasing_risk_2 <- tibble::tibble(x = 7.5, y = rep(5, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(4.5, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x, x = x - 0.13, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)


  # Create table with distribution of high risk values
  high_risk <- tibble::tibble(x = 7.5, y = rep(5.5, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(5, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x, x = x - 0.13, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  # Create table with distribution of high risk values
  ultra_risk <- tibble::tibble(x = 7.5, y = rep(7, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(5.5, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x, x = x - 0.13, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)


  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = increasing_risk_1,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
      size = 2
    ) +
    ggplot2::scale_color_gradient2(
      low = yellow, mid = yellow, high = red,
      midpoint = 3.5
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_segment(
      data = increasing_risk_2,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
      size = 2
    ) +
    ggplot2::scale_color_gradient2(
      low = red, mid = red, high = purple,
      midpoint = 4.5
    ) +

    ggplot2::scale_x_continuous(
      limits = c(0, 10)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 7)
    ) +
    # Add polar grid
    ggplot2::coord_polar(start = pi) +
    # Remove axis labels, grid and ticks
    ggplot2::theme_void() +
    # Remove legend
    ggplot2::theme(legend.position = "none")


  # Reset colour scale for next segment to be filled with new colour gradient
  p <- p + ggnewscale::new_scale_color() +
    ggplot2::geom_segment(
      data = high_risk,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
      size = 2
    ) +
    ggplot2::scale_color_gradient2(
      low = purple, mid = darkpurple, high = darkpurple,
      midpoint = 5.5
    )

  # Reset colour scale for next segment to be filled with new colour gradient
  p <- p + ggnewscale::new_scale_color() +
    ggplot2::geom_segment(
      data = safe_space,
      ggplot2::aes(
        x = x, xend = xend, y = y - 0.01, yend = yend - 0.01, color = y
      ),
      size = 2
    ) +
    ggplot2::scale_color_gradient2(
      low = "white",
      mid = "white",
      high = green,
      midpoint = 2
    )  +
    geom_segment(
      aes(x = 7.556, xend = 7.52, y = 2, yend = 5.61),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      size = 1.1,
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_segment(
      data = safe_space_arrow,
      ggplot2::aes(
        x = x, xend = xend, y = y - 0.02, yend = yend, alpha = y
      ),
      color = "white",
      size = 2
    ) +
    ggplot2::scale_alpha_continuous(
      range = c(0.4, 0),
    ) +
    geom_segment(
      aes(x = 7.54, xend = 7.54, y = 2.96, yend = 3.03),
      size = 3,
      color = "white"
    ) +
    geom_segment(
      aes(x = 7.53, xend = 7.53, y = 3.96, yend = 4.03),
      size = 3,
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.665,
        y = 2.5,
        label = "Safe Operating"
      ),
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.785,
        y = 2.555,
        label = "Space"
      ),
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.62,
        y = 3.515,
        label = "Zone of"
      ),
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.71,
        y = 3.525,
        label = "Uncertainty"
      ),
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.6,
        y = 4.445,
        label = "High Risk"
      ),
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.67,
        y = 4.455,
        label = "Zone"
      ),
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.58,
        y = 5.3,
        label = "Control"
      ),
      fontface = "italic",
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.64,
        y = 5.31,
        label = "Variable"
      ),
      fontface = "italic",
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 7.25, xend = 7.59, y = 3, yend = 3),
      size = 1.25,
      color = darkgreen,
      alpha = 0.8,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.25,
        y = 2.5,
        label = "Planetary Boundary"
      ),
      color = darkgreen,
      fontface = "bold",
      # alpha = 0.7,  # Add transparency
      size = 4  # Increase the size of the labels
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 7.27, xend = 7.575, y = 4, yend = 4),
      size = 1,
      size = 0.6,
      color = orange,
      alpha = 0.4,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.3,
        y = 3.65,
        label = "High Risk Line"
      ),
      color = orange,
      # alpha = 0.7,  # Add transparency
      size = 4.2  # Increase the size of the labels
    )

    final <- ggtrace::with_ggtrace(
      x = p + ggplot2::theme(aspect.ratio = .52),
      method = ggplot2::Layout$render,
      trace_steps = 5L,
      trace_expr = quote({
        panels <- lapply(
          panels,
          grid::editGrob,
          vp = grid::viewport(
            yscale = c(0.45, 1),
            xscale = c(0.45, 1)))
      }),
      out = "g"
    )

}