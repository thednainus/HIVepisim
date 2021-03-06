% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.initialize.R
\name{init_status_mig}
\alias{init_status_mig}
\title{Disease Status Initialization Module for netsim}
\usage{
init_status_mig(dat)
}
\arguments{
\item{dat}{Master list object containing a \code{networkDynamic} object and
other initialization information passed from \code{\link{netsim}}.}
}
\description{
This function sets the initial disease status on the
             network given the specified initial conditions.
}
\details{
This internal function sets, either randomly or deterministically, the nodes
that are infected at the starting time of network simulations, \eqn{t_1}.
If the number to be initially infected is passed, this function may set the
initial number infected based on the number specified, either as a a set of
random draws from a binomial distribution or as the exact number specified.
In either case, the specific nodes infected are a random sample from the
network. In contrast, a set of specific nodes may be infected by passing the
vector to \code{\link{netsim}}.

This module sets the time of infection for those nodes set infected
at the starting time of network simulations, \eqn{t_1}. For vital
dynamics models, the infection time for those nodes is a random draw from an
exponential distribution with the rate parameter defined by the
\code{di.rate} argument. For models without vital dynamics, the infection
time is a random draw from a uniform distribution of integers with a minimum
of 1 and a maximum of the number of time steps in the model. In both cases,
to set the infection times to be in the past, these times are multiplied by
-1, and 2 is added to allow for possible infection times up until step 2,
when the disease simulation time loop starts.
}
\seealso{
This is an initialization module for \code{\link{netsim}}.
}
\keyword{internal}
\keyword{netMod}
