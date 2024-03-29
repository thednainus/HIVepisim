% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_utils.R
\name{set_departures}
\alias{set_departures}
\title{Save id and time of individuals that departure the network}
\usage{
set_departures(dat, departures, at)
}
\arguments{
\item{dat}{Master list object containing a \code{networkDynamic} object and
other initialization information passed from \code{\link[EpiModel]{netsim}}.}

\item{departures}{ids of departures}

\item{at}{Current time step.}

\item{prefix}{Text for prefix to use when saving filename.}
}
\value{
the Master list object of network models (\code{dat})
}
\description{
It aims to save the ID and time of infected individual that
 departure the network via natural or HIV related cause.
}
