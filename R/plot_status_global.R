#' Plot the global status of planetary boundaries
#'
#' Plot line plots with the PB status over time  for
#' a scenario LPJmL run and derived planetary boundary statuses
#'
#' @param file_name character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen)
#'
#' @param x list with global output from calc_status
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
#'  plot_status(file_name = "./my_boundary_status.png",
#'                 x = list("land system change" = lsc_status,
#'                                    "nitrogen" = nitrogen_status,
#'                                    "bluewater" = water_status),
#'  legend = FALSE
#'  bg_col = NA)
#' }
#'
#' @md
#' @export

plot_status_global <- function(
                        x,
                        file_name = NULL, 
                        colors = c("safe zone" = '#e3f5ec',
                                   "increasing risk" = '#fcf8d8',
                                   "high risk" = '#f7e4dd'),
                        all_in_one = FALSE,
                        ncol = 2,
                        normalize = "safe"
                        ) {

  if (length(colors) != 3) {
    stop(paste0("Length of color vector (", length(colors), ") incorrect. ",
                 "Three colors have to be provided"))
  }

  normalize <- match.arg(normalize, c("safe", "increasing risk"))

  if (all_in_one == TRUE) {
    for (i in seq_len(length(x))) {
      x[[i]] <- as_risk_level(x[[i]], type = "continuous", normalize = normalize)
    }
  }
  
  data_tibble <- tidyr::as_tibble(x)
  data_tibble$years <- as.numeric(names(x[[1]]))
  data_tibble <- tidyr::pivot_longer(data_tibble, !years, names_to = "pb",
    values_to = "values")

  if (length(x) == 1 || all_in_one == TRUE) {
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
            height = 7 * ceiling(length(x) / ncol),
            units = "cm",
            res = 600,
            pointsize = 7)
      },
      `pdf` = {
        pdf(file_name,
            width = 10 * n_col / 2.54,
            height = 7 * ceiling(length(x) / ncol) / 2.54,
            pointsize = 7)
      }, {
        stop("File extension ", dQuote(file_extension), " not supported.")
      }
    )
  }

  if (all_in_one) {

    #holocene value
    holo_all <- NULL
    for (i in seq_len(length(x))) {
      holo_all <- abind(holo_all, attr(x[[i]], "thresholds")$holocene)
    }
    holo <- min(holo_all)
    pl_b <- attr(x[[1]], "thresholds")$pb
    h_risk <- attr(x[[1]], "thresholds")$highrisk
    plot <- ggplot2::ggplot(data_tibble, ggplot2::aes(x = years, y = values,
                                                      col = pb)) +
         ggplot2::geom_rect(ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf,
                                ymax =  pl_b), alpha = 1,
                            fill = colors[["safe zone"]], col = NA) +
         ggplot2::geom_rect(ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = pl_b,
                   ymax = max(values) + max(values) * 0.05),
                   alpha = 1, fill = colors[["increasing risk"]], col = NA)
    if (normalize == "safe") {
      plot <- plot + 
        ggplot2::scale_y_continuous(limits = c(holo, c(max(data_tibble$values) +
                                      max(data_tibble$values) * 0.05)),
                        expand = c(0, 0), breaks = c(holo, pl_b),
                        labels = c("holocene", "pb"))
    } else if (normalize == "increasing risk") {
      plot <- plot + 
        ggplot2::scale_y_continuous(limits = c(holo, c(max(data_tibble$values) +
                                     max(data_tibble$values) * 0.05)),
                       expand = c(0, 0), breaks = c(holo, pl_b, h_risk),
                       labels = c("holocene", "pb", "highrisk")) +
        ggplot2::geom_rect(ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = h_risk,
                                ymax =  Inf), alpha = 1,
                            fill = colors[["high risk"]], col = NA)
    }
    plot <- plot + ggplot2::geom_hline(yintercept = pl_b, linetype = 2, col = "grey") +
         ggplot2::geom_line() +
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
    numeric(length(x))
  for (i in seq_len(length(x))) {
    pb_thresh[i] <-  as.numeric(attr(x[[i]], "thresholds")[["pb"]])
    highrisk_thresh[i] <- as.numeric(
                            attr(x[[i]], "thresholds")[["highrisk"]]
                          )
    holocene[i] <- as.numeric(
                     attr(x[[i]], "thresholds")[["holocene"]]
                   )
    ylabel[i] <-  attr(x[[i]], "control variable")
    names(ylabel)[i] <- names(x)[i]
  }
  df_bg <- data.frame(
             pb = unique(data_tibble$pb),
             pb_thresh = pb_thresh,
             highrisk_thresh = highrisk_thresh,
             holocene = holocene,
             years = data_tibble$years[1:length(unique(data_tibble$pb))],
             values = data_tibble$values[1:length(unique(data_tibble$pb))]
           )
  min_years <- min(data_tibble$years)
  max_years <- max(data_tibble$years)
  darkgreen <- rgb(84. / 255., 156. / 255., 3. / 255., 1) #
  green <- rgb(134. / 255., 189. / 255., 36. / 255., 1)
  lightred <- rgb(244. / 255., 157. / 255., 136. / 255., 1)
  red <- rgb(197. / 255., 18. / 255., 22. / 255., 1)
  yellow <- rgb(251. / 255., 206. / 255., 0. / 255., 1)
  purple <- rgb(117. / 255., 24. / 255., 81. / 255., 1)


  # plot background
  plot <- ggplot2::ggplot() +
    ggpattern::geom_rect_pattern(data = df_bg,
                                 ggplot2::aes(xmin = min_years,
                                              xmax = max_years,
                                              ymin = holocene,
                                              ymax = pb_thresh),
                                 pattern = "gradient",
                                 pattern_fill = green,
                                 pattern_fill2 = darkgreen,
                                 alpha = 0.8) +
    ggpattern::geom_rect_pattern(data = df_bg,
                                 ggplot2::aes(xmin = min_years,
                                              xmax = max_years,
                                              ymin = pb_thresh,
                                              ymax = highrisk_thresh),
                                 pattern = "gradient",
                                 pattern_fill = yellow,
                                 pattern_fill2 = red,
                                 alpha = 0.8)  +
    ggpattern::geom_rect_pattern(data = df_bg,
                                 ggplot2::aes(xmin = min_years
                                              xmax = max_years,
                                              ymin = highrisk_thresh,
                                              ymax = highrisk_thresh + (highrisk_thresh - pb_thresh) * 1.5), #nolint
                                 pattern = "gradient",
                                 p tern_fill = red,
                                 p tern_fill2 = purple,
                                 a ha = 0.8)  +
    ggpattern::geom_rect_pattern(data = df_bg,
                                 ggplot2::aes(xmin = min_years,
                                              xmax = max_years,
                                              ymin = highrisk_thresh + (highrisk_thresh - pb_thresh) * 1.5), #nolint
                                              ymax = Inf,
                                 pattern = "gradient",
                                 pattern_fill = purple,
                                 pattern_fill2 = purple,
                                 alpha = 0.8)
   #ggplot2::coord_cartesian(ylim = c(0, pb_thresh + 3 * (holocene - pb_thresh)))
  
  # plot pb status timeseries
  #TODO the Y axis shouls be flexible depending on the PB
  plot <- plot +
    ggplot2::geom_ribbon(data = data_tibble, ggplot2::aes(x = years,
                         ymin = values, ymax = Inf), fill = "white") +
     ggplot2::geom_hline(data = df_bg,
                         mapping = ggplot2::aes(yintercept = pb_thresh),
                         col = darkgreen) +
     ggplot2::geom_hline(data = df_bg,
                         mapping = ggplot2::aes(yintercept = highrisk_thresh),
                         col = lightred, linewidth = 0.25) +
     ggplot2::geom_line(data = data_tibble,
                        ggplot2::aes(x = years, y = values, group = pb),
                        col = "gray31") +
     ggplot2::facet_wrap(~ pb, ncol = ncol, scales = "free",
                         labeller = ggplot2::as_labeller(ylabel),
                         strip.position = "left") +
     ggplot2::geom_label(data = data_tibble,
                         mapping = ggplot2::aes(x = -Inf, y = Inf,
                                                label = pb),
                         hjust = 0.09, vjust = 0.8, size = 3.5,
                         label.padding = ggplot2::unit(0.75, "lines")
                         ) +
     ggplot2::theme_classic(base_line_size = 0.25, base_rect_size = 0.25) +
     ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
     ggplot2::ylab("") +
     ggplot2::theme(strip.background = ggplot2::element_blank(),
           strip.placement = "outside") +
     ggplot2::theme(panel.border = ggplot2::element_rect(colour = "#6b6767",
                                                fill = NA, size = 0.5)) +
     ggplot2::scale_y_continuous(expand = c(0, NA))

  }
  print(plot)
  if (!is.null(file_name)) dev.off()
}
