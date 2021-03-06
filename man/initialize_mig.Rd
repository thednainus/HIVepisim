% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.initialize.R
\name{initialize_mig}
\alias{initialize_mig}
\title{Initialization: netsim Module}
\usage{
initialize_mig(x, param, init, control, s)
}
\arguments{
\item{x}{An \code{EpiModel} object of class \code{\link{netest}}.}

\item{param}{An \code{EpiModel} object of class \code{\link{param.net}}.}

\item{init}{An \code{EpiModel} object of class \code{\link{init.net}}.}

\item{control}{An \code{EpiModel} object of class \code{\link{control.net}}.}

\item{s}{Simulation number, used for restarting dependent simulations.}
}
\description{
This function initializes the master \code{dat} object on which
             data are stored, simulates the initial state of the network, and
             simulates disease status and other attributes.
}
\keyword{internal}
