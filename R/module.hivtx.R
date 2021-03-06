
#' @title Treatment Module
#'
#' @description Module function for anti-retroviral treatment (ART) initiation and
#'              adherence over time.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' Persons enter into the simulation with one of three ART "patterns": never
#' tested, tested but never treated, and treated. This module initiates ART
#' for treatment naive persons in the latter type, and then cycles them on
#' and off treatment. ART initiation, non-adherence, and restarting are all
#' stochastically simulated based on binomial statistical models.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{tx.status},
#' \code{tx.init.time}, \code{cuml.time.on.tx}, \code{cuml.time.off.tx} attributes.
#' Summary statistics for the average cumulative time on and off treatment
#' are calculated and stored on \code{dat$epi}.
#'
#'
#' @export
#'
hivtx_msm <- function(dat, at) {

  # Attributes ----
  status <- get_attr(dat, "status")
  tx.status <- get_attr(dat, "tx.status")
  diag.status <- get_attr(dat, "diag.status")
  cuml.time.on.tx <- get_attr(dat, "cuml.time.on.tx")
  cuml.time.off.tx <- get_attr(dat, "cuml.time.off.tx")
  tx.period.first <- get_attr(dat, "tx.period.first")
  tx.period.last <- get_attr(dat, "tx.period.last")
  tx.init.time <- get_attr(dat, "tx.init.time")

  # Parameters -----
  tx.init.prob <- get_param(dat, "tx.init.prob")
  tx.halt.prob <- get_param(dat, "tx.halt.prob")
  tx.reinit.prob <- get_param(dat, "tx.reinit.prob")

  art_start <- get_param(dat, "art_start")


  if (at >= art_start){
    ## Initiation -----
    tx.init.elig <- which(status == "i"  &
                           tx.status == 0  &
                           diag.status == 1  &
                           cuml.time.on.tx == 0 )
    rates <- tx.init.prob
    tx.init <- tx.init.elig[rbinom(length(tx.init.elig), 1, rates) == 1]

    ## Halting
    tx.halt.elig <- which(tx.status == 1 & cuml.time.on.tx > 0)
    rates.halt <- tx.halt.prob
    tx.halt <- tx.halt.elig[rbinom(length(tx.halt.elig), 1, rates.halt) == 1]

    ## Restarting
    tx.reinit.elig <- which(tx.status == 0 & cuml.time.on.tx > 0)
    rates.reinit <- tx.reinit.prob
    tx.reinit <- tx.reinit.elig[rbinom(length(tx.reinit.elig), 1, rates.reinit) == 1]

    ## Update Attributes
    #browser
    tx.status[tx.init] <- 1
    tx.status[tx.halt] <- 0
    tx.status[tx.reinit] <- 1

    # Save ART init ----
    if (dat$control$save.stats == TRUE){
      if (length(tx.init) > 0) {
        #browser()
        dat <- set_art_init(dat, at, tx.init)
      }

      if (length(tx.halt) > 0){
        dat <- set_art_halt(dat, at, tx.halt)
      }

      if (length(tx.reinit) > 0){
        dat <- set_art_reinit(dat, at, tx.reinit)
      }
    }

    cuml.time.on.tx[which(tx.status == 1)] <- cuml.time.on.tx[which(tx.status == 1)] + 1
    cuml.time.off.tx[which(tx.status == 0)] <- cuml.time.off.tx[which(tx.status == 0)] + 1

    tx.init.time[tx.init] <- at

    idsSetPeriod <- union(tx.init, tx.reinit)
    tx.period.first[idsSetPeriod] <- at
    tx.period.last[idsSetPeriod] <- at

    idsContPeriod <- setdiff(which(tx.status == 1), idsSetPeriod)
    tx.period.last[idsContPeriod] <- at

    # Output ----
    dat <- set_attr(dat, "tx.status", tx.status)
    dat <- set_attr(dat, "cuml.time.on.tx", cuml.time.on.tx)
    dat <- set_attr(dat, "cuml.time.off.tx", cuml.time.off.tx)
    dat <- set_attr(dat, "tx.period.first", tx.period.first)
    dat <- set_attr(dat, "tx.period.last", tx.period.last)
    dat <- set_attr(dat, "tx.init.time", tx.init.time)

    dat <- set_epi(dat, "mean.tx.on", at, mean(cuml.time.on.tx, na.rm = TRUE))
    dat <- set_epi(dat, "mean.tx.off", at, mean(cuml.time.off.tx, na.rm = TRUE))
  } else{
    dat <- set_epi(dat, "mean.tx.on", at, 0)
    dat <- set_epi(dat, "mean.tx.off", at, 0)
  }



  return(dat)
}
