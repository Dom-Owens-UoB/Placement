% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{fredmd}
\alias{fredmd}
\title{Get FRED-MD dataset}
\usage{
fredmd(
  file = "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv",
  date_start = NULL,
  date_end = NULL,
  transform = TRUE
)
}
\arguments{
\item{file}{string address to FREDMD data}

\item{date_start}{Date denoting start of period, of format YYYY-MM-01}

\item{date_end}{Date denoting end of period, of format YYYY-MM-01}

\item{transform}{Boolean, apply stationarity transforms}
}
\value{
data.frame of FRED-MD data
}
\description{
Get FRED-MD dataset
}
\examples{
fred_data <- fredmd()
}
