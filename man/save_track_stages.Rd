% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_utils.R
\name{save_track_stages}
\alias{save_track_stages}
\title{Save time and IDs of nodes and their stages of HIV infection.}
\usage{
save_track_stages(dat, prefix = NULL)
}
\arguments{
\item{dat}{Master list object containing a \code{networkDynamic} object and
other initialization information passed from \code{\link[EpiModel]{netsim}}.}

\item{prefix}{Text for prefix to use when saving filename.}
}
\value{
it writes a csv file at a specific location at the computer
}
\description{
Whenever an individual change the stage of HIV infection, this
  function will save the information on track_stages.
}
\details{
This function will save information to get the stage of infection at time
of sampling.


If a prefix is not provided, csv file will be saved as stages.csv
}
