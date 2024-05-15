#' Plot the global status of planetary boundaries
#'
#' Plot line plots with the PB status over time for
#' a scenario LPJmL run and derived planetary boundary statuses. Legend can be
#' plotted seperately based on the plot_legend() function
#'
#' @param x list with global output from calc_status
#'
#' @param filename character string providing file name (including directory
#'        and file extension). Defaults to NULL (plotting to screen and
#'        returning plot object for further customization)
#'
#' @param all_in_one boolean, if TRUE, all PB stati will be normalized and
#'        plotted in one panel
#'
#' @param ncol number of plot columns (only relevant if more than one pb is
#'        plotted and all_in_one = FALSE)
#'
#' @param normalize see [`as_risk_level()`] for details. Default set to
#'        "increasing risk". Only relevant if all_in_one = TRUE
#'
#' @param high_risk numeric, specify the normalized PB status value for the
#' upper end of the risk color scale (default 3.5)
#'
#' @examples
#' \dontrun{
#' plot_status_global(
#'   filename = "./my_boundary_status.png",
#'   x = status_output,
#'   all_in_one = FALSE,
#'   ncol = 2
#' )
#' }
#'
#' @md
#' @export
plot_status_global <- function(
    x,
    filename = NULL,
    all_in_one = FALSE,
    ncol = 2,
    normalize = "increasing risk",
    high_risk = 3.5) {

  normalize <- match.arg(normalize, c("safe", "increasing risk"))

  if (all_in_one == TRUE) {
    x <- as_risk_level(x, type = "continuous", normalize = normalize)
  }

  # please R CMD check for use of dplyr syntax
  long_name <- numeric(length(x))
  for (i in seq_len(length(x))) {
    class(x[[i]]) <- "numeric"
    long_name[i] <- attr(x[[i]], "long_name")
  }
  data_tibble <- tidyr::as_tibble(x)
  data_tibble$years <- as.numeric(names(x[[1]]))
  data_tibble <- tidyr::pivot_longer(
    data_tibble,
    !years,
    names_to = "pb",
    values_to = "values"
  )
  # define order of PBs
  data_tibble$pb <- factor(data_tibble$pb,
                           levels = c("lsc", "biosphere", "bluewater",
                                      "greenwater", "nitrogen"))


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
      ggplot2::geom_rect(
        fill = green,
        ggplot2::aes(
          xmin = -Inf, xmax = max(data_tibble$years),
          ymin = -Inf, ymax = pl_b
        ), alpha = 0.9
      )

    # please R CMD check for use of ggplot syntax
    if (normalize == "safe") {
      plot <- plot +
        ggplot2::scale_y_continuous(
          limits = c(holo, c(max(data_tibble$values) +
                               max(data_tibble$values) * 0.05)),
          expand = c(0, 0), breaks = c(holo, pl_b),
          labels = c("Holocene", "PB")
        ) +
        ggpattern::geom_rect_pattern(
          ggplot2::aes(
            xmin = -Inf, xmax = Inf,
            ymin = pl_b,
            ymax = max(data_tibble$values) +
              max(data_tibble$values) * 0.05
          ),
          pattern = "gradient",
          pattern_fill = yellow,
          pattern_fill2 = yellow,
          pattern_alpha = 0.5,
          col = NA
        )
    } else if (normalize == "increasing risk") {

      # create dataframe for plotting of gradient background filling
      df_inc_risk <- tibble::tibble()
      factor <- high_risk - 2 # 2 = upper end of increasing risk scale
      temp <- tibble::tibble(
        vals = seq(pl_b, (h_risk + (h_risk - pl_b) * factor),
          length.out = 500 # adjusted to not affect color scale
        )
      ) %>%
        dplyr::mutate(
          xend = max(data_tibble$years),
          x = min(data_tibble$years),
          yend = vals,
          col = seq(1, 0, length.out = length(vals))
        ) %>%
        dplyr::rename(y = vals)
      df_inc_risk <- dplyr::bind_rows(df_inc_risk, temp)

      # plot background
      plot <- plot + ggplot2::geom_segment(
        data = df_inc_risk,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = col),
        size = 0.5,
        show.legend = FALSE
      ) +
        ggplot2::scale_colour_viridis_c(
          option = "inferno",
          direction = 1,
          begin = 0.05,
          end = 0.85
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, c(max(data_tibble$values) +
                            max(data_tibble$values) * 0.05)),
          oob = scales::squish,
          expand = c(0, NA),
          breaks = c(holo, pl_b, h_risk),
          labels = c("Holocene", "PB", "High risk")
        )
    }
    # add horizontal lines for thresholds
    plot <- plot +
      ggplot2::geom_hline(
        yintercept = pl_b,
        linetype = 1,
        col = darkgreen,
        linewidth = 0.45
      )
    if (normalize == "increasing risk") {
      plot <- plot +
        ggplot2::geom_hline(
          yintercept = h_risk,
          linetype = 1,
          col = lightred,
          linewidth = 0.45
        )
    }
    # add trajectories for different pbs
    plot <- plot +
      ggnewscale::new_scale_colour() +
      # draw a slightly thicker black line to distinguish trajectories from background
      ggplot2::geom_line(
        data = data_tibble,
        ggplot2::aes(x = years, y = values, group = pb),
        size = 0.68, color = "black", alpha = 0.8
      ) +
      ggplot2::geom_line(
        data = data_tibble,
        ggplot2::aes(x = years, y = values, col = pb),
        size = 0.5,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_grey(start = 0.3, end = 0.7) +
      # add labels with the long names of the PBs
      # This will force the correct position of the link's right end
      ggrepel::geom_text_repel(
        data = data_tibble %>% dplyr::filter(years == max(years)),
        ggplot2::aes(x = years,
                     y = values,
                     label = gsub("^.*$", " ", long_name)),
        segment.square = TRUE,
        segment.color = "black", # "#6b6767",
        segment.alpha = 0.8,
        segment.size = 0.4,
        box.padding = 0.3,
        point.padding = 0.1,
        nudge_x = 0,
        nudge_y = -0.2, # change vertical alignment
        force = 0.5,
        hjust = 0.2,
        direction = "y",
        min.segment.length = 0.1,
        na.rm = TRUE,
        xlim = c(max(data_tibble$years) + 3,
                 max(data_tibble$years) + 11),
        ylim = c(0, max(data_tibble$values))
      ) +
      ggrepel::geom_text_repel(
        data = data_tibble %>% dplyr::filter(years == max(years)),
        ggplot2::aes(x = years,
                     y = values,
                     label = paste0("  ", long_name)),
        segment.alpha = 0, # This will 'hide' the link
        segment.square = TRUE,
        box.padding = 0.3,
        point.padding = 0.1,
        nudge_x = 0,
        nudge_y = -0.2, # change vertical alignment
        force = 0.5,
        hjust = 0.2,
        direction = "y",
        na.rm = TRUE,
        xlim = c(max(data_tibble$years), max(data_tibble$years) + 90),
        ylim = c(0, max(data_tibble$values)),
        min.segment.length = 0.1,
        size = 3,
        alpha = 0.8
      ) +
      ggplot2::theme_classic(base_line_size = 0.25, base_rect_size = 0.25) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        legend.position = "none"
      ) + # nolint:object_usage_linter
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          colour = "#6b6767",
          fill = NA, size = 0.5
        )
      ) +
      ggplot2::theme(axis.title = ggplot2::element_blank()) +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(, 3.2, , , "cm")
      ) +
      ggplot2::coord_cartesian(
        xlim = c(min(data_tibble$years), max(data_tibble$years)),
        expand = FALSE,
        clip = "off"
      )
  } else {
    # all_in_one = FALSE
    # create dataframe for plotting of pb specific background filling
    # (needed for safe zone and beyond the increasing risk value)
    pb_thresh <- highrisk_thresh <- holocene <- ylabel <- max_y <- long_name <-
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
      ylabel[i] <- paste0(
        attr(x[[i]], "control_variable"),
        " (", attr(x[[i]], "unit"), ")"
      )
      names(ylabel)[i] <- names(x)[i]
      long_name[i] <- paste0("  ", attr(x[[i]], "long_name"))
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

    # create dataframe for pb-specific background filling in the increasing risk
    # zone
    df_inc_risk <- tibble::tibble()
    factor <- high_risk - 2 # 2 = upper end of increasing risk scale
    for (i in seq_len(length(x))) {
      temp <- tibble::tibble(
        vals = seq(
          (highrisk_thresh[i] +
            (highrisk_thresh[i] - pb_thresh[i]) * factor), pb_thresh[i],
          length.out = 500
        )
      ) %>%
        dplyr::mutate(
          xend = min_years,
          x = max_years,
          yend = vals,
          pb = unique(data_tibble$pb)[i],
          col = seq(1, 0, length.out = length(vals))
        ) %>%
        dplyr::rename(y = vals)
      df_inc_risk <- dplyr::bind_rows(df_inc_risk, temp)
    }

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
        pattern_fill = green, # "#96c99e",
        pattern_fill2 = green, # "#96c99e",
        alpha = 0.8
      ) +
      ggplot2::geom_segment(
        data = df_inc_risk,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = col),
        size = 0.5,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_viridis_c(
        option = "inferno",
        direction = -1,
        begin = 0.15,
        end = 0.95
      ) +
      ggpattern::geom_rect_pattern(
        data = df_bg,
        ggplot2::aes(
          xmin = min_years,
          xmax = max_years,
          ymin = highrisk_thresh + (highrisk_thresh - pb_thresh) * factor
        ), # nolint
        ymax = Inf,
        pattern = "gradient",
        pattern_fill = "#000004FF",
        pattern_fill2 = "#000004FF",
        alpha = 0.8
      )

    # create scales for each pb
    upper_lim <- ifelse(max_y > highrisk_thresh, max_y, highrisk_thresh)
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
        col = lightred, linewidth = 0.45
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
        data = df_bg,
        mapping = ggplot2::aes(
          x = -Inf, y = Inf,
          label = long_name
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
      ggplot2::scale_y_continuous(expand = c(0, NA)) +
      ggplot2::scale_x_continuous(
        limits = range(data_tibble$years),
        expand = c(0, 0)
      )
  }

  if (length(x) == 1 || all_in_one == TRUE) {
    n_col <- 1
  } else {
    n_col <- ncol
  }

  if (is.null(filename)) {
    print(plot)
    return(plot)
  } else {
    width <- 9 * n_col
    if (all_in_one == FALSE) {
      height <- 7 * ceiling(length(x) / ncol)
    } else {
      height <- 6 * n_col * 0.75
    }
    ggplot2::ggsave(
      filename,
      plot,
      width = width,
      height = height,
      dpi = 600,
      units = "cm",
      pointsize = 7
    )
  }
}
