% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_utils.R
\name{save_stage}
\alias{save_stage}
\title{Save stage of HIV infection of nodes in the network at final step}
\usage{
save_stage(dat, prefix = NULL)
}
\arguments{
\item{dat}{Master list object containing a \code{networkDynamic} object and
other initialization information passed from \code{\link[EpiModel]{netsim}}.}

\item{prefix}{Text for prefix to use when saving filename.}
}
\value{
It writes a csv file to the specific location in the computer
}
\description{
At the final step of network simulation, it will save the IDs of
   infected nodes and their stage of HIV infection.
}
\details{
If a prefix is not provided, csv file will be saved as stage_and_IDs.csv
}
