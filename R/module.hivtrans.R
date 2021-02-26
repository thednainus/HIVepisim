#' @title Transmission Module
#'
#' @description Stochastically simulates disease transmission given the current
#'              state of the discordand edgelist. This function should be used
#'              when there is migration.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' This is the final substantive function that occurs within the time loop at
#' each time step. This function takes the discordant edgelist and calculates a
#' transmission probability for each row between dyads on the
#' network. After transmission events, individual-level attributes for the infected
#' persons are updated and summary statistics for incidence calculated.
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall are calculated and stored on \code{dat$epi}.
#'
#'
#' @export
#'
hivtrans_mig <- function(dat, at) {

  # Variables -----------------------------------------------------------

  # Attributes
  active <- get_attr(dat, "active")
  stage <- get_attr(dat, "stage")
  stage.time <- get_attr(dat, "stage.time")
  status <- get_attr(dat, "status")
  diag.status <- get_attr(dat, "diag.status")
  tx.status <- get_attr(dat, "tx.status")
  risk.group <- get_attr(dat, "risk.group")
  cuml.time.on.tx <- get_attr(dat, "cuml.time.on.tx")
  cuml.time.off.tx <- get_attr(dat, "cuml.time.off.tx")
  infTime <- get_attr(dat, "infTime")
  migrant <- get_attr(dat, "migrant")
  origin <- get_attr(dat, "origin")

  # Parameters ------
  # baseline for transmission rate
  trans.r <- get_param(dat, "trans.r")
  time.unit <- get_param(dat, "time.unit")

  # transmission risk ratio by stage of HIV infection
  ws0 <- get_param(dat, "ws0")
  ws1 <- get_param(dat, "ws1")
  ws2 <- get_param(dat, "ws2")
  ws3 <- get_param(dat, "ws3")
  ws4 <- get_param(dat, "ws4")

  # Transmission risk ratio
  # by care status (undiagnosed; diagnosed and untreated; and diagnosed an treated)
  wc1 <- get_param(dat, "wc1")
  wc2 <- get_param(dat, "wc2")
  wc3 <- get_param(dat, "wc3")

  # Transmission risk ratio
  # by risk group
  wr1 <- get_param(dat, "wr1")
  wr2 <- get_param(dat, "wr2")


  ## Find infected nodes ##
  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)


  ## Initialize default incidence at 0 ##
  nInf <- 0
  nInf.pop1 <- 0
  nInf.pop2 <- 0

  if (nElig > 0 && nElig < nActive) {

    ## Look up discordant edgelist ##
    del <- discord_edgelist(dat, at, network = 1)


    ## If any discordant pairs, proceed ##
    if (!(is.null(del))) {

      # Infection duration to at
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      del$status <- status[del$inf]
      del$risk.group <- risk.group[del$inf]
      del$stage <- stage[del$inf]
      del$tx.status <- tx.status[del$inf]
      del$cuml.time.on.tx <- cuml.time.on.tx[del$inf]
      del$cuml.time.off.tx <- cuml.time.off.tx[del$inf]
      del$diag.status <- diag.status[del$inf]
      del$infOrigin <- origin[del$inf]
      del$infMigrant <- migrant[del$inf]

      del$susOrigin <- origin[del$sus]
      del$susMigrant <- migrant[del$sus]
      del$susStatus <- status[del$sus]


      # Set parameters on discordant edgelist data frame
      del$trans.r <- rep(trans.r, length(del$inf))

      # stages of HIV infection (following Cori et al. 2015)
      del$trans.r[del$stage == 0] <- del$trans.r[del$stage == 0] * ws0
      del$trans.r[del$stage == 1] <- del$trans.r[del$stage == 1] * ws1
      del$trans.r[del$stage == 2] <- del$trans.r[del$stage == 2] * ws2
      del$trans.r[del$stage == 3] <- del$trans.r[del$stage == 3] * ws3
      del$trans.r[del$stage == 4] <- del$trans.r[del$stage == 4] * ws4

      # Treatment status and tested status
      #is not tested
      del$trans.r[del$diag.status == 0] <- del$trans.r[del$diag.status == 0] * wc1
      # is off treatment
      del$trans.r[del$tx.status == 0 & del$cuml.time.off.tx > 0] <-
        del$trans.r[del$tx.status == 0 & del$cuml.time.off.tx > 0] * wc2
      #is tested but not on treatment
      del$trans.r[del$diag.status == 1 & del$tx.status == 0 & del$cuml.time.on.tx == 0] <-
        del$trans.r[del$diag.status == 1 & del$tx.status == 0 & del$cuml.time.on.tx == 0] * wc2
      #is on treatment
      del$trans.r[del$tx.status == 1 & del$cuml.time.on.tx > 0] <-
        del$trans.r[del$tx.status == 1 & del$cuml.time.on.tx > 0] * wc3

      # risk groups
      del$trans.r[del$risk.group == 1] <- del$trans.r[del$risk.group == 1] * wr1
      del$trans.r[del$risk.group == 2] <- del$trans.r[del$risk.group == 2] * wr2


      # Transmission from infected person --------------------------------------

      transmit <- rbinom(nrow(del), 1, del$trans.r)
      del <- del[which(transmit == 1), ]

      # Look up new ids if any transmissions occurred
      ## Update Nodal Attr
      idsNewInf <- unique(del$sus)
      status <- get_attr(dat, "status")
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)

      nInf <- length(idsNewInf)
      nInf.pop1 <- sum(origin[idsNewInf] == "region")
      nInf.pop2 <- sum(origin[idsNewInf] == "global")
    }# end some discordant edges condition
  } # end some active discordant nodes condition

  if (nInf > 0) {

    # Attributes of newly infected
    status[idsNewInf] <- "i"
    dat <- set_attr(dat, "status", status)
    infTime[idsNewInf] <- at
    dat <- set_attr(dat, "infTime", infTime)
    stage[idsNewInf] <- 0
    dat <- set_attr(dat, "stage", stage)
    stage.time[idsNewInf] <- 0
    dat <- set_attr(dat, "stage.time", stage.time)
    diag.status[idsNewInf] <- 0
    dat <- set_attr(dat, "diag.status", diag.status)
    tx.status[idsNewInf] <- 0
    dat <- set_attr(dat, "tx.status", tx.status)
    cuml.time.on.tx[idsNewInf] <- 0
    dat <- set_attr(dat, "cuml.time.on.tx", cuml.time.on.tx)
    cuml.time.off.tx[idsNewInf] <- 0
    dat <- set_attr(dat, "cuml.time.off.tx", cuml.time.off.tx)


     # Attributes of transmitter
    transmitter <- as.numeric(del$inf)
    tab.trans <- table(transmitter)
    uni.trans <- as.numeric(names(tab.trans))
    count.trans <- get_attr(dat, "count.trans")
    count.trans[uni.trans] <- count.trans[uni.trans] + as.numeric(tab.trans)
    dat <- set_attr(dat, "count.trans", count.trans)
  }


  #Output ----

  # Save transmission matrix
  if (dat$control$save.transmat == TRUE){
    if (nInf > 0) {
      dat <- set_transmat(dat, del, at)
    }
  }

  ## Save incidence vector
  dat <- set_epi(dat, "incid.all", at, nInf)
  dat <- set_epi(dat, "incid.pop1", at, nInf.pop1)
  dat <- set_epi(dat, "incid.pop2", at, nInf.pop2)

  return(dat)
}
