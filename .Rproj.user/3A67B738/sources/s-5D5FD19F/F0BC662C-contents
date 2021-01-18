#' @title Arrivals Module
#'
#' @description Module function for arrivals into the sexually active
#'              population. Arrivals are modelled independetly for population 1
#'              (which will have origin = "region" and related to the MSM population in
#'              San Diego) and population 2 (which will have origin = global
#'              representing the global population of MSM).
#'
#' @inheritParams EpiModelHIV::aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions. For each new
#' entry, a set of attributes is added for that node, and the nodes are added onto
#' the network objects. Only attributes that are a part of the network model
#' formulae are updated as vertex attributes on the network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm
#' @export
#'
arrivals_mig <- function(dat, at){

  arrival.age <- get_param(dat, "arrival.age")

  # Variables ---------------------------------------------------------------
  # Arrival rate for population 1 (San Diego population)
  a1.rate <- get_param(dat, "a1.rate")
  index1 <- at - 1
  nOld1 <- get_epi(dat, "num.pop1", index1)
  nArrivals1 <- 0

  # Arrival rate for population 2 (Global population)
  a2.rate <- get_param(dat, "a2.rate")
  index2 <- at - 1
  nOld2 <- get_epi(dat, "num.pop2", index2)
  nArrivals2 <- 0

  # Add Nodes ---------------------------------------------------------------
  # For population 1 (San Diego population)
  if (nOld1 > 0) {
    nArrivals1 <- rbinom(1, nOld1, a1.rate)
    if (nArrivals1 > 0) {
      #dat <- append_core_attr(dat, at, nArrivals1)
      dat <- append_core_attr_mig(dat, at, nArrivals1)
      dat <- append_attr(dat, "status", "s", nArrivals1)
      dat <- append_attr(dat, "infTime", NA, nArrivals1)
      dat <- append_attr(dat, "age", arrival.age, nArrivals1)
      dat <- append_attr(dat, "diag.status", 0, nArrivals1)
      dat <- append_attr(dat, "count.trans", 0, nArrivals1)
      risk.group <- sample(1:2, nArrivals1, replace = TRUE)
      dat <- append_attr(dat, "risk.group", risk.group, nArrivals1)
      dat <- append_attr(dat, "num.neg.tests", 0, nArrivals1)
      dat <- append_attr(dat, "origin", "region", nArrivals1)
      dat <- append_attr(dat, "migrant", 1, nArrivals1)
      dat <- append_attr(dat, "stage", NA, nArrivals1)
      dat <- append_attr(dat, "stage.time", NA, nArrivals1)
      dat <- append_attr(dat, "aids.time", NA, nArrivals1)
      dat <- append_attr(dat, "diag.stage", NA, nArrivals1)
      dat <- append_attr(dat, "diag.time", NA, nArrivals1)
      dat <- append_attr(dat, "last.neg.test", NA, nArrivals1)
      dat <- append_attr(dat, "tx.status", NA, nArrivals1)
      dat <- append_attr(dat, "cuml.time.on.tx", NA, nArrivals1)
      dat <- append_attr(dat, "cuml.time.off.tx", NA, nArrivals1)
      dat <- append_attr(dat, "tx.period.first", NA, nArrivals1)
      dat <- append_attr(dat, "tx.period.last", NA, nArrivals1)
      dat <- append_attr(dat, "tx.init.time", NA, nArrivals1)
      #dat <- append_attr(dat, "migrationTime", NA, nArrivals1)
    }
  }


  # For population 2 (San Diego population)
  if (nOld2 > 0) {
    nArrivals2 <- rbinom(1, nOld2, a2.rate)
    if (nArrivals2 > 0) {
      #dat <- append_core_attr(dat, at, nArrivals2)
      dat <- append_core_attr_mig(dat, at, nArrivals2)
      dat <- append_attr(dat, "status", "s", nArrivals2)
      dat <- append_attr(dat, "infTime", NA, nArrivals2)
      dat <- append_attr(dat, "age", arrival.age, nArrivals2)
      dat <- append_attr(dat, "diag.status", 0, nArrivals2)
      dat <- append_attr(dat, "count.trans", 0, nArrivals2)
      risk.group <- sample(1:2, nArrivals2, replace = TRUE)
      dat <- append_attr(dat, "risk.group", risk.group, nArrivals2)
      dat <- append_attr(dat, "num.neg.tests", 0, nArrivals2)
      dat <- append_attr(dat, "origin", "global", nArrivals2)
      dat <- append_attr(dat, "migrant", 2, nArrivals2)
      dat <- append_attr(dat, "stage", NA, nArrivals2)
      dat <- append_attr(dat, "stage.time", NA, nArrivals2)
      dat <- append_attr(dat, "aids.time", NA, nArrivals2)
      dat <- append_attr(dat, "diag.stage", NA, nArrivals2)
      dat <- append_attr(dat, "diag.time", NA, nArrivals2)
      dat <- append_attr(dat, "last.neg.test", NA, nArrivals2)
      dat <- append_attr(dat, "tx.status", NA, nArrivals2)
      dat <- append_attr(dat, "cuml.time.on.tx", NA, nArrivals2)
      dat <- append_attr(dat, "cuml.time.off.tx", NA, nArrivals2)
      dat <- append_attr(dat, "tx.period.first", NA, nArrivals2)
      dat <- append_attr(dat, "tx.period.last", NA, nArrivals2)
      dat <- append_attr(dat, "tx.init.time", NA, nArrivals2)
      #dat <- append_attr(dat, "migrationTime", NA, nArrivals2)
    }
  }



  # Output ------------------------------------------------------------------
  dat <- set_epi(dat, "a1.flow", at, nArrivals1)
  dat <- set_epi(dat, "a2.flow", at, nArrivals2)

  return(dat)
}
