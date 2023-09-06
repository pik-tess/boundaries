#' Plot the global status of planetary boundaries
#'
#' Plot line plots with the PB status over time  for
#' a scenario LPJmL run and derived planetary boundary statuses
#'
#' @param file_name character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen)
#'
#' @param status_data list with global output from calc_status
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
#' @examples
#' \dontrun{
#'  plot_status(file_name = "./my_boundary_status.png",
#'                 status_data = list("land system change" = lsc_status,
#'                                    "nitrogen" = nitrogen_status,
#'                                    "bluewater" = water_status),
#'  legend = FALSE
#'  bg_col = NA)
#' }
#'
#' @md
#' @export

plot_status_global <- function(file_name = NULL,
                        status_data = NULL,
                        colors = c("safe zone" = '#e3f5ec',
                                   "increasing risk" = '#fcf8d8',
                                   "high risk" = '#f7e4dd'),
                        all_in_one = FALSE,
                        ncol = 2
                        ) {

  if (length(colors) != 3) {
    stop(paste0("Length of color vector (", length(colors), ") incorrect. ",
                 "Three colors have to be provided"))
  }
  if (all_in_one == TRUE) {
    for (i in seq_len(length(status_data))) {
      status_data[[i]] <- as_risk_level(status_data[[i]], type = "continuous")
    }
  }

  data_tibble <- tidyr::as_tibble(status_data)
  data_tibble$years <- as.numeric(names(status_data[[1]]))
  data_tibble <- tidyr::pivot_longer(data_tibble, !years, names_to = "pb",
    values_to = "values")



  if (length(status_data) == 1 || all_in_one == TRUE) {
    n_col <- 1
  } else {
    n_col <- ncol
  }
  if (!is.null(file_name)) {
    file_extension <- strsplit(file_name, split = "\\.")[[1]][-1] %>%
                      tail(1)
    switch(file_extension,
      `png` = {
        png(file_name,
            width = 9 * n_col,
            height = 7 * ceiling(length(status_data) / ncol),
            units = "cm",
            res = 600,
            pointsize = 7)
      },
      `pdf` = {
        pdf(file_name,
            width = 9 * n_col / 2.54,
            height = 7 * ceiling(length(status_data) / ncol) / 2.54,
            pointsize = 7)
      }, {
        stop("File extension ", dQuote(file_extension), " not supported.")
      }
    )
  }

  if (all_in_one == TRUE) {
    plot <- ggplot2::ggplot(data_tibble, ggplot2::aes(x = years, y = values,
                                                      col = pb)) +
         ggplot2::geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf,
                                ymax = 1), alpha = 1,
                            fill = colors[["safe zone"]], col = NA) +
         ggplot2::geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1,
                   ymax = max(values) + max(values) * 0.05),
                   alpha = 1, fill = colors[["increasing risk"]], col = NA) +
         ggplot2::geom_hline(yintercept = 1, linetype = 2, col = "grey") +
         ggplot2::geom_line() +
         ggplot2::scale_y_continuous(limits = c(0, c(max(data_tibble$values) +
                                          max(data_tibble$values) * 0.05)),
                            expand = c(0, 0), breaks = c(0, 1),
                            labels = c("holocene", "pb")) +
         ggplot2::scale_x_continuous(limits = range(data_tibble$years),
                            expand = c(0, 0)) +
         ggplot2::theme_classic(base_line_size = 0.25, base_rect_size = 0.25) +
         ggplot2::theme(legend.title = ggplot2::element_blank()) +
         ggplot2::theme(panel.border = ggplot2::element_rect(colour = "#6b6767",
                                                    fill = NA, size = 0.5)) +
         ggplot2::theme(axis.title = element_blank())
  } else {
  # create dataframe for plotting of pb specific background filling
    pb_thresh <- highrisk_thresh <- holocene <- ylabel <-
      numeric(length(status_data))
    for (i in seq_len(length(status_data))) {
      pb_thresh[i] <-  as.numeric(attr(status_data[[i]], "thresholds")[["pb"]])
      highrisk_thresh[i] <- as.numeric(
                              attr(status_data[[i]], "thresholds")[["highrisk"]]
                            )
      holocene[i] <- as.numeric(
                       attr(status_data[[i]], "thresholds")[["holocene"]]
                     )
      ylabel[i] <-  attr(status_data[[i]], "control variable")
      names(ylabel)[i] <- names(status_data)[i]
    }

    df_bg <- data.frame(
               pb = unique(data_tibble$pb),
               pb_thresh = pb_thresh,
               highrisk_thresh = highrisk_thresh,
               holocene = holocene,
               years = data_tibble$years[1:length(unique(data_tibble$pb))],
               values = data_tibble$values[1:length(unique(data_tibble$pb))]
             )
    # plot background
    plot <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = df_bg, ggplot2::aes(
                                        xmin = -Inf, xmax = +Inf,
                                        ymin = holocene, ymax = pb_thresh),
                         fill = "#e3f5ec", alpha = 0.8) +
      ggplot2::geom_rect(data = df_bg, ggplot2::aes(xmin = -Inf, xmax = +Inf,
                                        ymin = pb_thresh,
                                        ymax = highrisk_thresh),
                         fill = "#fcf8d8", alpha = 0.8) +
      ggplot2::geom_rect(data = df_bg, ggplot2::aes(
                                            xmin = -Inf, xmax = +Inf,
                                            ymin = highrisk_thresh,
                                            ymax = +Inf),
                           fill = "#f7e4dd", alpha = 0.8) +
      ggplot2::geom_hline(data = df_bg,
                          mapping = ggplot2::aes(yintercept = pb_thresh),
                          linetype = 2, col = "grey") +
      ggplot2::scale_y_continuous(expand = c(0, NA))
     #ggplot2::coord_cartesian(ylim = c(0, pb_thresh + 3 * (holocene - pb_thresh)))
    
    # plot pb status timeseries
    #TODO y axis name needs to be pb specific (as attribute of status_data)
    #TODO PB itself should be at the same position in all plots
    plot <- plot +
         ggplot2::geom_line(data = data_tibble,
                            ggplot2::aes(x = years, y = values, group = pb),
                            col = "#8e8a8a") +
         ggplot2::facet_wrap(~ pb, ncol = ncol, scales = "free",
                             labeller = as_labeller(ylabel),
                             strip.position = "left") +
         ggplot2::geom_label(data = data_tibble,
                             mapping = ggplot2::aes(x = -Inf, y = Inf,
                                                    label = pb),
                             hjust = 0.09, vjust = 0.8, size = 3.5,
                             label.padding = unit(0.75, "lines")
                             ) +
         ggplot2::theme_classic(base_line_size = 0.25, base_rect_size = 0.25) +
         ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
         ggplot2::ylab("") +
         theme(strip.background = element_blank(),
               strip.placement = "outside") +
         ggplot2::theme(panel.border = ggplot2::element_rect(colour = "#6b6767",
                                                    fill = NA, size = 0.5))

  }
  print(plot)
  if (!is.null(file_name)) dev.off()
}
