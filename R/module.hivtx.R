
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
  origin <- get_attr(dat, "origin")
  active <- get_attr(dat, "active")
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
    ## Initiation: "region" population -----
    tx.init.elig_pop1 <- which(active == 1 &
                                 status == "i" &
                                 tx.status == 0 &
                                 diag.status == 1 &
                                 cuml.time.on.tx == 0 &
                                 origin == "region")
    tx.init_pop1 <- tx.init.elig_pop1[rbinom(length(tx.init.elig_pop1), 1,
                                             tx.init.prob) == 1]

    ## Halting: "region" population
    tx.halt.elig_pop1 <- which(active == 1 & tx.status == 1 & cuml.time.on.tx > 0 &
                                 origin == "region")
    tx.halt_pop1 <- tx.halt.elig_pop1[rbinom(length(tx.halt.elig_pop1), 1, tx.halt.prob) == 1]

    ## Restarting: "region" population
    tx.reinit.elig_pop1 <- which(active == 1 & tx.status == 0 & cuml.time.on.tx > 0 &
                                   origin == "region")
    tx.reinit_pop1 <- tx.reinit.elig_pop1[rbinom(length(tx.reinit.elig_pop1), 1,
                                                 tx.reinit.prob) == 1]

    ## Initiation: "global" population -----
    tx.init.elig_pop2 <- which(active == 1 &
                                 status == "i" &
                                 tx.status == 0 &
                                 diag.status == 1 &
                                 cuml.time.on.tx == 0 &
                                 origin == "global")
    tx.init_pop2 <- tx.init.elig_pop2[rbinom(length(tx.init.elig_pop2), 1,
                                             tx.init.prob) == 1]

    ## Halting: "global" population
    tx.halt.elig_pop2 <- which(active == 1 & tx.status == 1 & cuml.time.on.tx > 0 &
                                 origin == "global")
    tx.halt_pop2 <- tx.halt.elig_pop2[rbinom(length(tx.halt.elig_pop2), 1, tx.halt.prob) == 1]

    ## Restarting: "global" population
    tx.reinit.elig_pop2 <- which(active == 1 & tx.status == 0 & cuml.time.on.tx > 0 &
                                   origin == "global")
    tx.reinit_pop2 <- tx.reinit.elig_pop2[rbinom(length(tx.reinit.elig_pop2), 1,
                                                 tx.reinit.prob) == 1]


    ## Update Attributes
    #browser
    if (length(tx.init_pop1) > 0){
      tx.status[tx.init_pop1] <- 1
    }

    if(length(tx.halt_pop1) > 0) {
      tx.status[tx.halt_pop1] <- 0
    }

    if (length(tx.reinit_pop1) > 0) {
      tx.status[tx.reinit_pop1] <- 1
    }





    if (length(tx.init_pop2) > 0){
      tx.status[tx.init_pop2] <- 1
    }

    if(length(tx.halt.elig_pop2) > 0) {
      tx.status[tx.halt.elig_pop2] <- 0
    }

    if (length(tx.reinit.elig_pop2) > 0) {
      tx.status[tx.reinit.elig_pop2] <- 1
    }



    # Save ART init ----
    if (dat$control$save.stats == TRUE){
      if (length(tx.init_pop1) > 0 | length(tx.init_pop2) > 0) {
        #browser()
        dat <- set_art_init(dat, at, c(tx.init_pop1, tx.init_pop2))
      }

      if (length(tx.halt_pop1) > 0 | length(tx.halt_pop2) > 0){
        dat <- set_art_halt(dat, at, c(tx.halt_pop1, tx.halt_pop2))
      }

      if (length(tx.reinit_pop2) > 0 | length(tx.reinit_pop1) > 0){
        dat <- set_art_reinit(dat, at, c(tx.halt_pop1_pop1, tx.halt_pop1_pop2))
      }
    }

    cuml.time.on.tx[which(tx.status == 1)] <- cuml.time.on.tx[which(tx.status == 1)] + 1
    cuml.time.off.tx[which(tx.status == 0)] <- cuml.time.off.tx[which(tx.status == 0)] + 1

    tx.init.time[c(tx.init_pop1, tx.init_pop2)] <- at

    #browser()
    idsSetPeriod <- c(tx.init_pop1, tx.init_pop2, tx.reinit_pop1, tx.reinit_pop2)
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
