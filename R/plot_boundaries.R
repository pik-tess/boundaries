# install.packages("remotes")
# remotes::install_github("yjunechoe/ggtrace")

library(ggplot2)
library(grid)
library(tidyr)
library(ggnewscale)
library(ggtrace)

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

data_path <- (
  "/p/projects/open/Johanna/boundaries/R/r_out/r_data/global_timeseries_porkka_50.RData" # nolint
)

data_path_ma <- (
  "/p/projects/open/Johanna/boundaries/R/r_out/r_data/global_timeseries_avg_10_porkka.RData"
)


data <- load(data_path_ma) |>
  get()

for (i in seq_len(length(data))) {
  data[[i]] <- boundaries::as_risk_level(
    data[[i]],
    type = "continuous",
    normalize = "increasing risk"
  )
}

data <- lapply(data, function(x){
  x[x < 0 | is.nan(x) | is.na(x)] <- 0
  x[x > 3] <- 3
  x
})

# Create demo data for now
final_data <- tibble(
  x = seq_along(data),
  y = rep(3, length(data)),
  start_year = min(as.integer(names(data[[names(data)[1]]]))),
  end_year = max(as.integer(names(data[[names(data)[1]]]))),
  mid_year = round((start_year + end_year) / 2),
  # years = as.integer(names(data[[names(data)[1]]])),
  name = names(data),
  status = data
)

# Create ggproto class to plot time series of boundaries into wedge
TSBoundaryStatus <- ggplot2::ggproto( # nolint
  "TSBoundaryStatus",
  Stat,
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
  Stat,
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
  layer(
    data = data,
    mapping = mapping,
    stat = if (negative) TSBoundaryStatusNegative else TSBoundaryStatus, # nolint
    geom = GeomPolygon,
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
safe_space <- tibble(x = 1:5, y = rep(1, 5)) %>%
  dplyr::mutate(vals = purrr::map(y, ~seq(0, .x, by = 0.01))) %>%
  tidyr::unnest(cols = c(vals)) %>%
  dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
  dplyr::select(-y) %>%
  dplyr::rename(y = vals)

# Create table with distribution of uncertainty zone values
uncertainty_zone <- tibble(x = 1:5, y = rep(2, 5)) %>%
  dplyr::mutate(vals = purrr::map(y, ~seq(1, .x, by = 0.01))) %>%
  tidyr::unnest(cols = c(vals)) %>%
  dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
  dplyr::select(-y) %>%
  dplyr::rename(y = vals)

# Create table with distribution of high risk values
high_risk <- tibble(x = 1:5, y = rep(3, 5)) %>%
  dplyr::mutate(vals = purrr::map(y, ~seq(2, .x, by = 0.01))) %>%
  tidyr::unnest(cols = c(vals)) %>%
  dplyr::mutate(xend = x + 0.44, x = x - 0.44, yend = vals) %>%
  dplyr::select(-y) %>%
  dplyr::rename(y = vals)


# Create ggplot, first with the segments that include the colour gradient for
#   each level of risk (safe space, uncertainty zone, high risk)
p <- ggplot() +
  geom_segment(
    data = high_risk,
    aes(x = x, xend = xend, y = y, yend = yend, color = y),
    size = 2
  ) +
  scale_color_gradient2(
    low = lightred, mid = red, high = purple,
    midpoint = max(high_risk$y) / 2
  ) +
  # Reset colour scale for next segment to be filled with new colour gradient
  ggnewscale::new_scale_color() +
  geom_segment(
    data = uncertainty_zone,
    aes(x = x, xend = xend, y = y, yend = yend, color = y),
    size = 2
  ) +
  scale_color_gradient2(
    low = yellow, mid = yellow, high = orange,
    midpoint = max(uncertainty_zone$y) / 2
  ) +
  # Reset colour scale for next segment to be filled with new colour gradient
  ggnewscale::new_scale_color() +
  geom_segment(
    data = safe_space,
    aes(x = x, xend = xend, y = y, yend = yend, color = y),
    size = 2
  ) +
  scale_color_gradient2(
    low = lightgreen,
    mid = green,
    high = darkgreen,
    midpoint = max(safe_space$y) / 2
  ) +
  # Hide segment area outside of boundaries with white polygon
  #   (complementary part of wedge)
  geom_boundaries(
    data = final_data,
    aes(x = x, y = y, group = x, status = status),
    fill = "white",
    negative = TRUE
  ) +
  # Draw outline around actual boundaries time series wedge
  geom_boundaries(
    data = final_data,
    aes(x = x, y = y, group = x, status = status),
    color = "black",
    fill = "transparent"
  ) +
  # Limit x-axis to a half circle
  scale_x_continuous(
    limits = c(0, pi^2)
  ) +
  # Add polar grid
  coord_polar(start = pi * 1.392) +
  # Remove axis labels, grid and ticks
  theme_void() +
  # Remove legend
  theme(legend.position = "none")

p <- p +
  geom_segment(
    data = final_data,
    aes(x = x - 0.45, xend = x - 0.45, y = 0, yend = 3.12),
    size = 0.8,
    color = "black",
    alpha = 0.7,  # Add transparency
    linetype = 1
  ) +
  geom_segment(
    data = final_data,
    aes(x = x, xend = x, y = 3, yend = 3.12),
    size = 0.8,
    color = "black",
    alpha = 0.7,  # Add transparency
    linetype = 1
  ) +
  geom_segment(
    data = final_data,
    aes(x = x + 0.45, xend = x + 0.45, y = 0, yend = 3.12),
    size = 0.8,
    color = "black",
    alpha = 0.7,  # Add transparency
    linetype = 1
  )

# Add lines/segments to indicate the different levels of risk
p <- p +
  geom_segment(
    data = final_data,
    aes(x = x - 0.45, xend = x + 0.45, y = 1, yend = 1),
    size = 0.8,
    color = darkgreen,
    linetype = 1
  ) +
  geom_segment(
    data = final_data,
    aes(x = x - 0.45, xend = x + 0.45, y = 2, yend = 2),
    size = 0.8,
    color = orange,
    linetype = 1
  ) +
  geom_segment(
    data = final_data,
    aes(x = x - 0.45, xend = x + 0.45, y = 3, yend = 3),
    size = 0.8,
    color = purple,
    linetype = 1
  )

# Add lines/segments to seperate the different boundaries
p <- p +
  geom_segment(
    aes(x = 0.53, xend = 0.53, y = 0, yend = 3.5),
    size = 1.2,
    color = "black",
    linetype = 1
  ) +
  geom_segment(
    aes(x = 1.5, xend = 1.5, y = 0, yend = 3.5),
    size = 1.2,
    color = "black",
    linetype = 1
  ) +
  geom_segment(
    aes(x = 2.5, xend = 2.5, y = 0, yend = 3.5),
    size = 1.2,
    color = "black",
    linetype = 1
  ) +
  geom_segment(
    aes(x = 3.5, xend = 3.5, y = 0, yend = 3.5),
    size = 1.2,
    color = "black",
    linetype = 1
  ) +
  geom_segment(
    aes(x = 4.5, xend = 4.5, y = 0, yend = 3.5),
    size = 1.2,
    color = "black",
    linetype = 1
  ) +
  geom_segment(
    aes(x = 5.47, xend = 5.47, y = 0, yend = 3.5),
    size = 1.2,
    color = "black",
    linetype = 1
  )

# Add texts on the plot
p <- p +
  geom_text(
    data = final_data,
    aes(x = x - 0.43, y = 3.25, label = start_year),
    color = "black",
    alpha = 0.7,  # Add transparency
    size = 4  # Increase the size of the labels
  ) +
  geom_text(
    data = final_data,
    aes(x = x, y = 3.25, label = mid_year),
    color = "black",
    alpha = 0.7,  # Add transparency
    size = 4  # Increase the size of the labels
  ) +
  geom_text(
    data = final_data,
    aes(x = x + 0.43, y = 3.25, label = end_year),
    color = "black",
    alpha = 0.7,  # Add transparency
    size = 4  # Increase the size of the labels
  )

# Add texts on the plot
p <- p +
  geom_label(
    data = final_data,
    aes(x = x, y = 3.8, label = name),
    color = "black",
    size = 6  # Increase the size of the labels
  )

# Reformat ggplot for better half circle display
final <- ggtrace::with_ggtrace(
  x = p + theme(aspect.ratio = .52),
  method = Layout$render,
  trace_steps = 5L,
  trace_expr = quote({
    panels <- lapply(panels, editGrob, vp = viewport(yscale = c(0.48, 1)))
  }),
  out = "g"
)
 
ggsave("test.png", final, width = 16, height = 9, dpi = 300)
