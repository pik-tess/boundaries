% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_risk_level.R
\name{as_risk_level}
\alias{as_risk_level}
\title{Convert status of control variable to risk level}
\usage{
as_risk_level(control_variable, type = "continuous", normalize = "safe")
}
\arguments{
\item{control_variable}{output array from calc_* with the status of the
control variable, incl. pb thresholds as attribute}

\item{type}{character string to define whether to return risk level as
continuous (normalized so that 0 = holocene state and
1 = planetary boundary; >1 = transgressed) or discrete variable
(0 = no PB status assessed, 1 = safe, 2 = increasing risk,
3 = high risk)}

\item{normalize}{character string to define normalization, either "safe"
(normalized from holocene to pb = the safe zone) or
"increasing risk" (normalized from pb to high risk level =
increasing risk zone if the pb status is > pb, otherwise normalized
from holocene to pb). Only used if type set to "continuous"}
}
\description{
Convert status of control variable to planetary boundary risk level
(safe, increasing risk, high risk), based on the output from calc_*
}
\examples{
\dontrun{
as_risk_level(
  control_variable = biosphere_status,
  type = "discrete"
)
}

}
