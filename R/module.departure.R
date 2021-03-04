#' @title Depature Module
#'
#' @description Module function for simulating both general and disease-related
#'              departures, including deaths, among population members.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and AIDS-mortality rate.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#'
#' @keywords module msm
#' @export
#'
departure_mig <- function(dat, at) {

  ## General departures
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  origin <- get_attr(dat, "origin")
  exitTime <- get_attr(dat, "exitTime")
  rates.all <- get_param(dat, "asmr")
  rates.aids <- get_param(dat, "aids.mr")

  idsDpt.pop1 <- NULL
  idsDpt.aids.pop1 <- NULL

  idsDpt.pop2 <- NULL
  idsDpt.aids.pop2 <- NULL

  # Departures (not HIV-related: population 1) --------------------------------------------------
  nDepartures.pop1 <- 0
  idsElig.pop1 <- which(active == 1 & origin == "region")
  nElig.pop1 <- length(idsElig.pop1)
  if (nElig.pop1 > 0) {
    vecDepartures.pop1 <- which(rbinom(nElig.pop1, 1, rates.all) == 1)
    if (length(vecDepartures.pop1) > 0) {
      idsDpt.pop1 <- idsElig.pop1[vecDepartures.pop1]
      nDepartures.pop1 <- length(idsDpt.pop1)
      active[idsDpt.pop1] <- 0
      exitTime[idsDpt.pop1] <- at
    }
  }

  # Departures (not HIV-related: population 2) --------------------------------------------------
  nDepartures.pop2 <- 0
  idsElig.pop2 <- which(active == 1 & origin == "global")
  nElig.pop2 <- length(idsElig.pop2)
  if (nElig.pop2 > 0) {
    vecDepartures.pop2 <- which(rbinom(nElig.pop2, 1, rates.all) == 1)
    if (length(vecDepartures.pop2) > 0) {
      idsDpt.pop2 <- idsElig.pop2[vecDepartures.pop2]
      nDepartures.pop2 <- length(idsDpt.pop2)
      active[idsDpt.pop2] <- 0
      exitTime[idsDpt.pop2] <- at
    }
  }


  # AIDS-related departures (pop1) -----------------------------------------------------
  nDepartures.aids.pop1 <- 0
  idsElig.aids.pop1 <- which(active == 1 & stage == 4 & origin == "region")
  nElig.aids.pop1 <- length(idsElig.aids.pop1)
  if (nElig.aids.pop1 > 0) {
    vecDepartures.aids.pop1 <- which(rbinom(nElig.aids.pop1, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop1) > 0) {
      idsDpt.aids.pop1 <- idsElig.aids.pop1[vecDepartures.aids.pop1]
      nDepartures.aids.pop1 <- length(idsDpt.aids.pop1)
      active[idsDpt.aids.pop1] <- 0
      exitTime[idsDpt.aids.pop1] <- at
    }
  }

  # AIDS-related departures (pop2) -----------------------------------------------------
  nDepartures.aids.pop2 <- 0
  idsElig.aids.pop2 <- which(active == 1 & stage == 4 & origin == "global")
  nElig.aids.pop2 <- length(idsElig.aids.pop2)
  if (nElig.aids.pop2 > 0) {
    vecDepartures.aids.pop2 <- which(rbinom(nElig.aids.pop2, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop2) > 0) {
      idsDpt.aids.pop2 <- idsElig.aids.pop2[vecDepartures.aids.pop2]
      nDepartures.aids.pop2 <- length(idsDpt.aids.pop2)
      active[idsDpt.aids.pop2] <- 0
      exitTime[idsDpt.aids.pop2] <- at
    }
  }


  # counting the number of departures by natural causes
  # that was HIV positive
  idsDepAll.pop1 <- unique(c(idsDpt.pop1, idsDpt.aids.pop1))
  depHIV.pop1 <- intersect(idsDepAll.pop1, which(status == "i"))
  ndepHIV.pop1 <- length(depHIV.pop1)


  idsDepAll.pop2 <- unique(c(idsDpt.pop2, idsDpt.aids.pop2))
  depHIV.pop2 <- intersect(idsDepAll.pop2, which(status == "i"))
  ndepHIV.pop2 <- length(depHIV.pop2)


  #Cumulative R0 calculations
  if (at == 2) {
    dat$temp$R0_pop1 <- NA
    dat$temp$R0_pop2 <- NA
  }
  if (length(depHIV.pop1) > 0) {
    newR0_pop1 <- dat$attr$count.trans[depHIV.pop1]
    dat$temp$R0_pop1 <- c(dat$temp$R0_pop1, newR0_pop1)
  }
  if (length(depHIV.pop2) > 0) {
    newR0_pop2 <- dat$attr$count.trans[depHIV.pop2]
    dat$temp$R0_pop2 <- c(dat$temp$R0_pop2, newR0_pop2)
  }


  # Output ------------------------------------------------------------------

  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  dat <- set_epi(dat, "dall_pop1.flow", at, nDepartures.pop1)
  dat <- set_epi(dat, "daids_pop1.flow", at, nDepartures.aids.pop1)
  dat <- set_epi(dat, "dhiv_pop1.flow", at, ndepHIV.pop1)

  dat <- set_epi(dat, "dall_pop2.flow", at, nDepartures.pop2)
  dat <- set_epi(dat, "daids_pop2.flow", at, nDepartures.aids.pop2)
  dat <- set_epi(dat, "dhiv_pop2.flow", at, ndepHIV.pop2)

  return(dat)
}

