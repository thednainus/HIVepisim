#' @title Prevalence Calculations within Time Steps for 2 populations
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules (arrival and departure)
#' as they depend on vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#'
#' @export
#'
prevalence_mig <- function(dat, at) {

  #Attributes ----

  active <- get_attr(dat, "active")
  origin <- get_attr(dat, "origin")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  diag.status <- get_attr(dat, "diag.status")
  diag.stage <- get_attr(dat, "diag.stage")
  diag.time <- get_attr(dat, "diag.time")
  #diag.time <- abs(diag.time)
  aids.time <- get_attr(dat, "aids.time")
  infTime <- get_attr(dat, "infTime")
  age <- get_attr(dat, "age")
  tx.init.time <- get_attr(dat, "tx.init.time")
  last.neg.test <- get_attr(dat, "last.neg.test")

  if (at == 1) {
    dat$epi <- list()
    # set initial conditions for arrivals and migrations in the epi vector
    dat <- set_epi(dat, "a1.flow", at, 0)
    dat <- set_epi(dat, "a2.flow", at, 0)

    dat <- set_epi(dat, "nArrivals_mig1", at, 0)
    dat <- set_epi(dat, "nArrivals_mig2", at, 0)

    dat <- set_epi(dat, "dall.flow", at, 0)
    dat <- set_epi(dat, "daids.flow", at, 0)
    dat <- set_epi(dat, "dhiv.flow", at, 0)

    dat <- set_epi(dat, "mean.tx.on", at, 0)
    dat <- set_epi(dat, "mean.tx.off", at, 0)
    dat <- set_epi(dat, "newDx", at, 0)
    dat <- set_epi(dat, "new.aids.tot", at, 0)

    dat <- set_epi(dat, "incid.all", at, 0)
    dat <- set_epi(dat, "incid.pop1", at, 0)
    dat <- set_epi(dat, "incid.pop2", at, 0)
    dat <- set_epi(dat, "tot.tests", at, 0)
    dat <- set_epi(dat, "tot.neg.tests", at, 0)

  }

  # Pop Size / Demog ----
  # Total
  dat <- set_epi(dat, "num", at, sum(active == 1))
  dat <- set_epi(dat, "num.pop1", at, sum(active == 1 & origin == "region"))
  dat <- set_epi(dat, "num.pop2", at, sum(active == 1 & origin == "global"))



  #susceptibles
  dat <- set_epi(dat, "s.num.all", at, sum(active == 1 & status == "s", na.rm = TRUE))
  dat <- set_epi(dat, "s.num.pop1", at, sum(active == 1 & status == "s" & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "s.num.pop2", at, sum(active == 1 & status == "s" & origin == "global", na.rm = TRUE))


  #age
  dat <- set_epi(dat, "age.all", at, mean(age[active == 1], na.rm = TRUE))
  dat <- set_epi(dat, "age.pop1", at, mean(age[active == 1 & origin == "region"], na.rm = TRUE))
  dat <- set_epi(dat, "age.pop2", at, mean(age[active == 1 & origin == "global"], na.rm = TRUE))


  #infected
  dat <- set_epi(dat, "i.num.all", at, sum(active == 1 & status == "i", na.rm = TRUE))
  dat <- set_epi(dat, "i.num.pop1", at, sum(active == 1 & status == "i" & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "i.num.pop2", at, sum(active == 1 & status == "i" & origin == "global", na.rm = TRUE))



  # diagnostic
  dat <- set_epi(dat, "i.num.dx.all", at, sum(active == 1 & diag.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "i.num.dx.pop1", at, sum(active == 1 & diag.status == 1 & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "i.num.dx.pop2", at, sum(active == 1 & diag.status == 1 & origin == "global", na.rm = TRUE))


  # Prevalence
  dat <- set_epi(dat, "i.prev.all", at, dat$epi$i.num.all[at] / dat$epi$num[at])
  dat <- set_epi(dat, "i.prev.pop1", at, dat$epi$i.num.pop1[at] / dat$epi$num.pop1[at])
  dat <- set_epi(dat, "i.prev.pop2", at, dat$epi$i.num.pop2[at] / dat$epi$num.pop2[at])


  dat <- set_epi(dat, "i.prev.dx.all", at, dat$epi$i.num.dx.all[at] / dat$epi$num[at])
  dat <- set_epi(dat, "i.prev.dx.pop1", at, dat$epi$i.num.dx.pop1[at] / dat$epi$num.pop1[at])
  dat <- set_epi(dat, "i.prev.dx.pop2", at, dat$epi$i.num.dx.pop2[at] / dat$epi$num.pop2[at])


  # HIV stage
  dat <- set_epi(dat, "hstage0.all", at, sum(stage == 0 & diag.time >= 2 & active == 1, na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hstage0.pop1", at, sum(stage == 0 & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "hstage0.pop2", at, sum(stage == 0 & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE))



  dat <- set_epi(dat, "hstage1.all", at, sum(stage == 1 & diag.time >= 2 & active == 1, na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hstage1.pop1", at, sum(stage == 1 & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "hstage1.pop2", at, sum(stage == 1 & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE))


  dat <- set_epi(dat, "hstage2.all", at, sum(stage == 2 & diag.time >= 2 & active == 1, na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hstage2.pop1", at, sum(stage == 2 & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "hstage2.pop2", at, sum(stage == 2 & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE))



  dat <- set_epi(dat, "hstage3.all", at, sum(stage == 3 & diag.time >= 2 & active == 1, na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hstage3.pop1", at, sum(stage == 3 & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "hstage3.pop2", at, sum(stage == 3 & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE))

  dat <- set_epi(dat, "hstage.aids.all", at, sum(stage == 4 & diag.time >= 2 & active == 1, na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hstage.aids.pop1", at, sum(stage == 4 & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "region", na.rm = TRUE))
  dat <- set_epi(dat, "hstage.aids.pop2", at, sum(stage == 4 & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE) /
                   sum(status == "i" & diag.time >= 2 & active == 1 & origin == "global", na.rm = TRUE))


  return(dat)
}
