#' Plot the global status of planetary boundaries
#'
#' Plot line plots with the PB status over time  for
#' a scenario LPJmL run and derived planetary boundary statuses
#'
#' @param x list with global output from calc_status
#'
#' @param filename character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen)
#'
#' @param add_legend logical, specify whether a legend should be plotted
#'
#' @param colors named vector for definition of colors for plotting of
#'        the safe, increasingnrisk and high risk zone
#'
#' @param all_in_one boolean, if TRUE, all PB stati will be normalized to 0
#'         (holocene value) to 1 (pb) and plotted in one panel
#'
#' @param ncol number of plot columns (only relevant if more than one pb is
#'        plotted and all_in_one = FALSE)
#'
#' @param normalize character string to define normalization, either "safe"
#'        (normalized from holocene to pb = the safe zone) or
#'        "increasing risk" (normalized from pb to high risk level =
#'        increasing risk zone if the pb status is > pb, otherwise normalized
#'        from holocene to pb). Only used if all_on_one set to TRUE
#'
#' @examples
#' \dontrun{
#' plot_status(
#'   filename = "./my_boundary_status.png",
#'   x = list(
#'     "land system change" = lsc_status,
#'     "nitrogen" = nitrogen_status,
#'     "bluewater" = water_status
#'   ),
#'   legend = FALSE,
#'   bg_col = NA
#' )
#' }
#'
#' @md
#' @export
plot_status_global <- function(
    x,
    filename = NULL,
    add_legend = TRUE,
    colors = c(
      "bluewater" = "#377eb8", #' #4d4bb2', #"#d4d4d4", #
      "greenwater" = "#4daf4a", #' #41865c', "#b4b4b4", #
      "lsc" = "#ff7f00", #' #b77f06', "#909090",#
      "biosphere" = "#e41a1c", # '#b53a0d', "#636363", ##
      "nitrogen" = "#984ea3"
    ), #' #8a026a'), "#494848"), #
    all_in_one = FALSE,
    ncol = 2,
    normalize = "safe") {

  normalize <- match.arg(normalize, c("safe", "increasing risk"))

  if (all_in_one == TRUE) {
    x <- as_risk_level(x, type = "continuous", normalize = normalize)
  }

  # please R CMD check for use of dplyr syntax
  years <- NULL

  data_tibble <- tidyr::as_tibble(x)
  data_tibble$years <- as.numeric(names(x[[1]]))
  data_tibble <- tidyr::pivot_longer(
    data_tibble,
    !years,
    names_to = "pb",
    values_to = "values"
  )

  if (length(x) == 1 || all_in_one == TRUE) {
    n_col <- 1
  } else {
    n_col <- ncol
  }
  if (!is.null(filename)) {
    file_extension <- file_ext(filename)
    width <- 9 * n_col
    if (all_in_one == FALSE) {
      height <- 7 * ceiling(length(x) / ncol)
    } else {
      height <- 9 * n_col * 0.75
    }
    switch(file_extension,
      `png` = {
        grDevices::png(filename,
          width = width,
          height = height,
          units = "cm",
          res = 600,
          pointsize = 7
        )
      },
      `pdf` = {
        grDevices::pdf(filename,
          width = width / 2.54,
          height = height / 2.54,
          pointsize = 7
        )
      }
    )
  }

  if (all_in_one) {
    # holocene value
    holo_all <- NULL
    for (i in seq_len(length(x))) {
      holo_all <- abind::abind(holo_all, attr(x[[i]], "thresholds")$holocene)
    }
    holo <- min(holo_all)
    pl_b <- attr(x[[1]], "thresholds")$pb
    h_risk <- attr(x[[1]], "thresholds")$highrisk
    plot <- ggplot2::ggplot() +
      ggpattern::geom_rect_pattern(
        ggplot2::aes(
          xmin = -Inf, xmax = Inf,
          ymin = -Inf, ymax = pl_b
        ),
        pattern = "gradient",
        pattern_fill = "white",
        pattern_fill2 = green,
        pattern_alpha = 0.5,
        col = NA
      )

    # please R CMD check for use of ggplot syntax
    values <- pb <- NULL
    if (normalize == "safe") {
      plot <- plot +
        ggplot2::scale_y_continuous(
          limits = c(holo, c(max(data_tibble$values) +
                               max(data_tibble$values) * 0.05)),
          expand = c(0, 0), breaks = c(holo, pl_b),
          labels = c("holocene", "pb")
        ) +
        ggpattern::geom_rect_pattern(
          ggplot2::aes(
            xmin = -Inf, xmax = Inf,
            ymin = pl_b,
            ymax = max(values) +
              max(values) * 0.05
          ),
          pattern = "gradient",
          pattern_fill = yellow,
          pattern_fill2 = yellow,
          pattern_alpha = 0.5,
          col = NA
        )
    } else if (normalize == "increasing risk") {
      plot <- plot +
        ggplot2::coord_cartesian(
          ylim = c(holo, c(max(data_tibble$values) +
                             max(data_tibble$values) * 0.05)
          )
        ) +
        ggplot2::scale_y_continuous(
          expand = c(0, 0),
          breaks = c(holo, pl_b, h_risk),
          labels = c("holocene", "pb", "highrisk")
        ) +
        ggpattern::geom_rect_pattern(
          ggplot2::aes(
            xmin = -Inf, xmax = Inf,
            ymin = pl_b,
            ymax = h_risk
          ),
          pattern = "gradient",
          pattern_fill = yellow,
          pattern_fill2 = orange,
          pattern_alpha = 0.5,
          col = NA
        ) +
        ggpattern::geom_rect_pattern(
          ggplot2::aes(
            xmin = -Inf,
            xmax = Inf,
            ymin = h_risk,
            ymax =  (h_risk + (h_risk - pl_b) * 1.5)
          ),
          pattern = "gradient",
          pattern_fill = orange,
          pattern_fill2 = darkpurple,
          pattern_alpha = 0.5,
          col = NA
        )
    }

    plot <- plot +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
        fill = "white",
        col = NA,
        alpha = 0.68
      ) +
      ggplot2::geom_hline(
        yintercept = pl_b,
        linetype = 1,
        col = darkgreen
      )
    if (normalize == "increasing risk") {
      plot <- plot +
        ggplot2::geom_hline(
          yintercept = h_risk,
          linetype = 1,
          col = lightred,
          linewidth = 0.25
        )
    }
    plot <- plot +
      ggplot2::geom_line(
        data = data_tibble,
        mapping = ggplot2::aes(x = years, y = values, col = pb)
      ) +
      ggplot2::scale_color_manual(values = colors) +

      ggplot2::scale_x_continuous(
        limits = range(data_tibble$years),
        expand = c(0, 0)
      ) +
      ggplot2::theme_classic(base_line_size = 0.25, base_rect_size = 0.25) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) + # nolint:object_usage_linter
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          colour = "#6b6767",
          fill = NA, size = 0.5
        )
      ) +
      ggplot2::theme(axis.title = ggplot2::element_blank()) +
      ggplot2::theme(aspect.ratio = 1)
  } else {
    # create dataframe for plotting of pb specific background filling
    pb_thresh <- highrisk_thresh <- holocene <- ylabel <- max_y <-
      numeric(length(x))
    for (i in seq_len(length(x))) {
      pb_thresh[i] <- as.numeric(attr(x[[i]], "thresholds")[["pb"]])
      highrisk_thresh[i] <- as.numeric(
        attr(x[[i]], "thresholds")[["highrisk"]]
      )
      holocene[i] <- as.numeric(
        attr(x[[i]], "thresholds")[["holocene"]]
      )
      max_y[i] <- max(data_tibble$values[which(data_tibble$pb == names(x)[i])])
      ylabel[i] <- attr(x[[i]], "control variable")
      names(ylabel)[i] <- names(x)[i]
    }
    df_bg <- data.frame(
      pb = unique(data_tibble$pb),
      pb_thresh = pb_thresh,
      highrisk_thresh = highrisk_thresh,
      holocene = holocene,
      years = data_tibble$years[seq_along(unique(data_tibble$pb))],
      values = data_tibble$values[seq_along(unique(data_tibble$pb))]
    )
    min_years <- min(data_tibble$years)
    max_years <- max(data_tibble$years)

    # plot background
    plot <- ggplot2::ggplot() +
      ggpattern::geom_rect_pattern(
        data = df_bg,
        ggplot2::aes(
          xmin = min_years,
          xmax = max_years,
          ymin = holocene,
          ymax = pb_thresh
        ),
        pattern = "gradient",
        pattern_fill = "white",
        pattern_fill2 = green,
        alpha = 0.8
      ) +
      ggpattern::geom_rect_pattern(
        data = df_bg,
        ggplot2::aes(
          xmin = min_years,
          xmax = max_years,
          ymin = pb_thresh,
          ymax = highrisk_thresh
        ),
        pattern = "gradient",
        pattern_fill = yellow,
        pattern_fill2 = orange,
        alpha = 0.8
      ) +
      ggpattern::geom_rect_pattern(
        data = df_bg,
        ggplot2::aes(
          xmin = min_years,
          xmax = max_years,
          ymin = highrisk_thresh,
          ymax = highrisk_thresh + (highrisk_thresh - pb_thresh) * 1.5
        ), # nolint
        pattern = "gradient",
        pattern_fill = orange,
        pattern_fill2 = darkpurple,
        alpha = 0.8
      ) +
      ggpattern::geom_rect_pattern(
        data = df_bg,
        ggplot2::aes(
          xmin = min_years,
          xmax = max_years,
          ymin = highrisk_thresh + (highrisk_thresh - pb_thresh) * 1.5
        ), # nolint
        ymax = Inf,
        pattern = "gradient",
        pattern_fill = darkpurple,
        pattern_fill2 = darkpurple,
        alpha = 0.8
      )

    # create scales for each pb
    upper_lim <- ifelse(max_y > highrisk_thresh, max_y, highrisk_thresh)
    max(data_tibble$values) + max(data_tibble$values) * 0.05
    df_scales <- data.frame(
      pb =  unique(data_tibble$pb),
      ymin = holocene,
      ymax = upper_lim + 0.05 * upper_lim
    )
    df_scales <- split(df_scales, df_scales$pb)
    scales <- lapply(df_scales, function(x) {
      ggplot2::scale_y_continuous(
        limits = c(x$ymin, x$ymax),
        expand = c(0, 0)
      )
    })

    # plot pb status timeseries
    plot <- plot +
      ggplot2::geom_ribbon(data = data_tibble, ggplot2::aes(
        x = years,
        ymin = values, ymax = Inf
      ), fill = "white") +
      ggplot2::geom_hline(
        data = df_bg,
        mapping = ggplot2::aes(yintercept = pb_thresh),
        col = darkgreen
      ) +
      ggplot2::geom_hline(
        data = df_bg,
        mapping = ggplot2::aes(yintercept = highrisk_thresh),
        col = lightred, linewidth = 0.25
      ) +
      ggplot2::geom_line(
        data = data_tibble,
        ggplot2::aes(x = years, y = values, group = pb),
        col = "gray31"
      ) +
      ggplot2::facet_wrap(~pb,
        ncol = ncol, scales = "free",
        labeller = ggplot2::as_labeller(ylabel),
        strip.position = "left"
      ) +
      ggh4x::facetted_pos_scales(y = scales) +
      ggplot2::geom_label(
        data = data_tibble,
        mapping = ggplot2::aes(
          x = -Inf, y = Inf,
          label = pb
        ),
        hjust = 0.09, vjust = 0.8, size = 3.5,
        label.padding = ggplot2::unit(0.75, "lines")
      ) +
      ggplot2::theme_classic(base_line_size = 0.25, base_rect_size = 0.25) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
      ggplot2::ylab("") +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.placement = "outside"
      ) +
      ggplot2::theme(panel.border = ggplot2::element_rect(
        colour = "#6b6767",
        fill = NA, size = 0.5
      )) +
      ggplot2::scale_y_continuous(expand = c(0, NA))
  }
  print(plot)
  if (!is.null(filename)) grDevices::dev.off()
}
