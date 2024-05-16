#' Plot the legend for the normalized colors of PB statuses
#'
#' Plot a legend for the colors of PB statuses, normalized based on the size of
#' the increasing risk, for globally aggregated plots,
#' or spatially distributed maps
#'
#' @param filename character string providing file name (including directory
#' and file extension). Defaults to NULL (return
#' plot object for further adaptation)
#'
#' @param fontsize numeric specifying the size of the font to be used for
#' legend labels. Default set to 3.
#'
#' @examples
#' \dontrun{
#' plot_legend(
#'   filename = "./mylegend.png",
#' )
#' }
#'
#' @md
#' @export
plot_legend <- function(filename = NULL, fontsize = 3) {

  # please R CMD check for use of dplyr syntax
  vals <- y <- x <- NULL
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
  increasing_risk <- tibble::tibble(x = 7.5, y = rep(5.5, 5)) %>%
    dplyr::mutate(vals = purrr::map(y, ~seq(3, .x, by = 0.01))) %>%
    tidyr::unnest(cols = c(vals)) %>%
    dplyr::mutate(xend = x, x = x - 0.13, yend = vals) %>%
    dplyr::select(-y) %>%
    dplyr::rename(y = vals)

  # please R CMD check for use of dplyr syntax
  xend <- yend <- NULL
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = increasing_risk,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = y),
      size = 2
    ) +
    # Use inferno colour scale for colours beyond safe space
    ggplot2::scale_colour_viridis_c(
      option = "inferno",
      begin = 0.1,
      end = 0.95,
      direction = -1
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
      data = safe_space,
      ggplot2::aes(
        x = x, xend = xend, y = y - 0.01, yend = yend - 0.01, color = y
      ),
      size = 2
    ) +
    ggplot2::scale_color_gradient2(
      low = green,
      mid = green,
      high = green,
      midpoint = 2
    )  +
    ggplot2::geom_segment(
      data = safe_space,
      ggplot2::aes(x = x, xend = xend, y = 2.983, yend = 2.983),
      size = 1.5,
      color = "white",
      alpha = 0.1,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_segment(
      data = safe_space,
      ggplot2::aes(x = x, xend = xend, y = 2.987, yend = 2.987),
      size = 0.8,
      color = "white",
      alpha = 0.6,  # Add transparency
      linetype = 1
    ) +
    ggplot2::geom_segment(
      data = safe_space,
      ggplot2::aes(x = x, xend = xend, y = 2.9885, yend = 2.9885),
      size = 0.4,
      color = "white",
      alpha = 0.9,  # Add transparency
      linetype = 1
    ) +

    ggplot2::geom_segment(
      ggplot2::aes(x = 7.556, xend = 7.52, y = 1.975, yend = 5.61),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.25, "cm"), type = "closed"
      ),
      size = 1.1,
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 7.54, xend = 7.54, y = 2.96, yend = 3.03),
      size = 3,
      color = "white"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 7.53, xend = 7.53, y = 3.96, yend = 4.03),
      size = 3,
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.75,
        y = 2.5,
        label = "Safe Operating \nSpace"
      ),
      color = "black",
      lineheight = 1,
      # alpha = 0.7,  # Add transparency
      size = fontsize  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.68,
        y = 3.515,
        label = "Zone of \nUncertainty"
      ),
      lineheight = 1,
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = fontsize  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.64,
        y = 4.445,
        label = "High Risk \nZone",
        lineheight = 1
      ),
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = fontsize  # Increase the size of the labels
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 7.62,
        y = 5.3,
        label = "Control \nVariable"
      ),
      fontface = "italic",
      lineheight = 1,
      color = "black",
      # alpha = 0.7,  # Add transparency
      size = fontsize  # Increase the size of the labels
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
      size = fontsize  # Increase the size of the labels
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 7.27, xend = 7.575, y = 4, yend = 4),
      size = 1,
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
      size = fontsize  # Increase the size of the labels
    )

  final <- ggtrace::with_ggtrace(
    x = p + ggplot2::theme(aspect.ratio = .095),
    method = ggplot2::Layout$render,
    trace_steps = 5L,
    trace_expr = quote({
      panels <- lapply(
        panels,
        grid::editGrob,
        vp = grid::viewport(
          yscale = c(0.46, 0.54),
          xscale = c(0.55, 0.90)
        )
      )
    }),
    out = "g"
  )

  if (is.null(filename)) {
    return(final)
  } else {
    ggplot2::ggsave(filename, final, width = 10, height = 2)
  }
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
