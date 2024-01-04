#' Convert status of control variable to risk level
#'
#' Convert status of control variable to planetary boundary risk level
#' (safe, increasing risk, high risk), based on the output from calc_*
#'
#' @param control_variable output array from calc_* with the status of the
#'        control variable, incl. pb thresholds as attribute
#'
#' @param type character string to define whether to return risk level as
#'             continuous (normalized so that 0 = holocene state and
#'             1 = planetary boundary; >1 = transgressed) or discrete variable
#'             (0 = no PB status assessed, 1 = safe, 2 = increasing risk,
#'             3 = high risk)
#' 
#' @param normalize character string to define normalization, either "safe"
#'        (normalized from holocene to pb = the safe zone) or
#'        "increasing risk" (normalized from pb to high risk level =
#'        increasing risk zone if the pb status is > pb, otherwise normalized
#'        from holocene to pb). Only used if type set to "continuous"
#'
#' @examples
#' \dontrun{
#'  as_risk_level(control_variable = biosphere_status,
#'                type = "discrete")
#' }
#'
#' @md
#' @export

as_risk_level <- function(control_variable, type = "continuous", normalize = "safe") {

  type <- match.arg(type, c("continuous", "discrete"))

  thresholds <- attr(control_variable, "thresholds")

  if (type == "continuous") {
    if (normalize == "safe") {
      # holocene: 0, pb: 1
      # control variable status normalized to 0-1, 1 is the
      # threshold between safe and increasing risk, >1 is transgressed
      if (class(thresholds) == "list") {
        risk_level <- (control_variable - thresholds[["holocene"]]) /
                    (thresholds[["pb"]] - thresholds[["holocene"]])
        attr(risk_level, "thresholds") <-
          list(holocene = 0, pb = 1,
               highrisk = (thresholds[["highrisk"]] -
                           thresholds[["holocene"]]) /
                          (thresholds[["pb"]] -
                           thresholds[["holocene"]]))
        
      } else if (class(thresholds) == "array") {
        risk_level <- (control_variable - thresholds[[, , "holocene"]]) /
                    (thresholds[[, , "pb"]] - thresholds[[, , "holocene"]])
        attr(risk_level, "thresholds") <-
          list(holocene = 0, pb = 1,
               highrisk = (thresholds[[, , "highrisk"]] -
                           thresholds[[, , "holocene"]]) /
                          (thresholds[[, , "pb"]] -
                           thresholds[[, , "holocene"]]))
      }
    } else if (normalize == "increasing risk") {
      # holocene: 0, pb: 1, high risk level: 2
      # if the control variable status is > pb:  
      # control variable status normalized to 1-2, 1 is the
      # threshold between safe and increasing risk (pb), 2 is the threshold
      # between increasing risk and high risk zone
      # if control variable is < pb: 
      # normalized to 0 (holocene) to 1 (pb)
      if (class(thresholds) == "list") {
        risk_level <- control_variable
        risk_level[control_variable >= thresholds[["pb"]]] <-
          (control_variable[control_variable >= thresholds[["pb"]]] -
           thresholds[["pb"]]) /
          (thresholds[["highrisk"]] - thresholds[["pb"]]) + 1
        risk_level[control_variable < thresholds[["pb"]]] <-
          1 - (control_variable[control_variable < thresholds[["pb"]]] -
           thresholds[["pb"]]) /
          (thresholds[["holocene"]] - thresholds[["pb"]])
        attr(risk_level, "thresholds") <-
          list(holocene = 0, 
               pb = 1, highrisk = 2)
          #alternative, if no additional normalization from holocene to pb:
          #holocene = (thresholds[["holocene"]] -
          #                   thresholds[["pb"]]) / 
          #                  (thresholds[["highrisk"]] -
          #                   thresholds[["pb"]])
      } else if (class(thresholds) == "array") {
        risk_level <- control_variable
        risk_level[control_variable >= thresholds[[, , "pb"]]] <-
          (control_variable[control_variable >= thresholds[[, , "pb"]]] -
           thresholds[[, , "pb"]]) /
          (thresholds[[, , "highrisk"]] - thresholds[[, , "pb"]]) + 1
        risk_level[control_variable < thresholds[[, , "pb"]]] <-
          1 - (control_variable[control_variable < thresholds[[, , "pb"]]] -
           thresholds[[, , "pb"]]) /
          (thresholds[[, , "holocene"]] - thresholds[[, , "pb"]])
        attr(risk_level, "thresholds") <-
          list(holocene = 0, 
               pb = 1, highrisk = 2)
      }
    }
    

  } else if (type == "discrete") {
    # init array based on control_variable
    risk_level <- control_variable
    if (class(thresholds) == "list") {
      # high risk
      risk_level[control_variable >= thresholds[["highrisk"]]] <- 3
      # increasing risk
      risk_level[control_variable < thresholds[["highrisk"]] &
                 control_variable >= thresholds[["pb"]]] <- 2
      # safe zone
      risk_level[control_variable < thresholds[["pb"]]] <- 1
    } else if (class(thresholds) == "array") {
      # high risk
      risk_level[control_variable >= thresholds[, , "highrisk"]] <- 3
      # increasing risk
      risk_level[control_variable < thresholds[, , "highrisk"] &
                 control_variable >= thresholds[, , "pb"]] <- 2
      # safe zone
      risk_level[control_variable < thresholds[, , "pb"]] <- 1
      risk_level[is.na(thresholds[, , "pb"])] <- 0

    }
    # non applicable cells
    risk_level[is.na(control_variable)] <- 0
  }
  return(risk_level)
}