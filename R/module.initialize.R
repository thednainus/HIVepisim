
#' @title Initialization: netsim Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @export
#' @keywords internal
#'
initialize_mig <- function(x, param, init, control, s) {

  #browser()

  if (control$start == 1) {

    # Master Data List --------------------------------------------------------
    dat <- create_dat_object(param, init, control)
    dat$param <- param
    dat$init <- init
    dat$control <- control

    dat$nwparam <- list()
    dat$nwparam[[1]] <- x[-which(names(x) == "fit")]


    # Initial Network Simulation ----------------------------------------------

    if (get_control(dat, "resimulate.network") == TRUE) {
      nsteps <- 1
    } else {
      nsteps <- get_control(dat, "nsteps")
    }
    dat <- sim_nets_t1(x, dat, nsteps)


    # Network Parameters ------------------------------------------------------
    dat$nwparam <- list(x[-which(names(x) == "fit")])

    groups <- length(unique(get_vertex_attribute(dat$nw[[1]], "group")))
    dat <- set_param(dat, "groups", groups)


    # Nodal Attributes --------------------------------------------------------

    # Standard attributes
    #browser()
    num <- network.size(dat$nw[[1]])
    #dat <- append_core_attr(dat, 1, num)
    dat <- append_core_attr_mig(dat, 1, num)

    groups <- length(unique(get_vertex_attribute(dat$nw[[1]], "group")))
    dat <- set_param(dat, "groups", groups)

    ## Pull attr on nw to dat$attr
    dat <- copy_nwattr_to_datattr(dat)

    ## Store current proportions of attr
    nwterms <- get_network_term_attr(dat$nw[[1]])
    if (!is.null(nwterms)) {
      dat$temp$nwterms <- nwterms
      dat$temp$t1.tab <- get_attr_prop(dat, nwterms)
    }

    ## Infection Status and Time
    dat <- init_status_mig(dat)

    # Conversions for tergmLite
    tergmLite <- get_control(dat, "tergmLite")
    if (tergmLite == TRUE) {
      dat <- tergmLite::init_tergmLite(dat)
    }


    # Summary Stats -----------------------------------------------------------
    dat <- do.call(control[["prevalence.FUN"]], list(dat, at = 1))


    # Restart/Reinit Simulations ----------------------------------------------
  } else if (control$start > 1) {
    dat <- create_dat_object(param = x$param, control = control)

    dat$nw <- x$network[[s]]
    dat$nwparam <- x$nwparam
    if (is.null(dat$control$isTERGM)) {
      nwparam <- get_nwparam(dat)
      isTERGM <- all(nwparam$coef.diss$duration > 1)
      dat <- set_control(dat, "isTERGM", isTERGM)
    }
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$stats <- sapply(x$stats, function(var) var[[s]])
  }

  return(dat)
}


#' @title Disease Status Initialization Module for netsim
#'
#' @description This function sets the initial disease status on the
#'              network given the specified initial conditions.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and
#'        other initialization information passed from \code{\link{netsim}}.
#'
#' @details
#' This internal function sets, either randomly or deterministically, the nodes
#' that are infected at the starting time of network simulations, \eqn{t_1}.
#' If the number to be initially infected is passed, this function may set the
#' initial number infected based on the number specified, either as a a set of
#' random draws from a binomial distribution or as the exact number specified.
#' In either case, the specific nodes infected are a random sample from the
#' network. In contrast, a set of specific nodes may be infected by passing the
#' vector to \code{\link{netsim}}.
#'
#' This module sets the time of infection for those nodes set infected
#' at the starting time of network simulations, \eqn{t_1}. For vital
#' dynamics models, the infection time for those nodes is a random draw from an
#' exponential distribution with the rate parameter defined by the
#' \code{di.rate} argument. For models without vital dynamics, the infection
#' time is a random draw from a uniform distribution of integers with a minimum
#' of 1 and a maximum of the number of time steps in the model. In both cases,
#' to set the infection times to be in the past, these times are multiplied by
#' -1, and 2 is added to allow for possible infection times up until step 2,
#' when the disease simulation time loop starts.
#'
#' @seealso This is an initialization module for \code{\link{netsim}}.
#'
#' @export
#' @keywords netMod internal
#'
init_status_mig <- function(dat) {


  type <- get_control(dat, "type", override.null.error = TRUE)
  type <- if (is.null(type)) "None" else type

  nsteps <- get_control(dat, "nsteps")
  tergmLite <- get_control(dat, "tergmLite")
  groups <- get_param(dat, "groups")
  #browser()
  status.vector <- get_init(dat, "status.vector", override.null.error = TRUE)

  #depar.rates.all <- get_param(dat, "asmr")
  depar.rates.aids <- get_param(dat, "aids.mr")

  isTERGM <- get_control(dat, "isTERGM")

  # Variables ---------------------------------------------------------------

  #browser()
  i.num.all <- get_init(dat, "i.num.all", override.null.error = TRUE)

  num.all <- sum(get_attr(dat, "active") == 1)

  if (groups == 2) {
    group <- get_attr(dat, "group")
    if (!all(group %in% c(1, 2))) {
      stop(
        "When using the `group` attribute, the only authorized values",
        " are 1 and 2.\n",
        "The values found were: ", paste0(unique(group), collapse = ", ")
      )
    }

    i.num.g2 <- get_init(dat, "i.num.g2")
    if (type  == "SIR" && is.null(status.vector)) {
      r.num.g2 <- get_init(dat, "r.num.g2", override.null.error = TRUE)
    }
  } else {
    group <- rep(1, num.all)
  }

  statOnNw <- "status" %in% dat$temp$nwterms

  # Status ------------------------------------------------------------------

  ## Status passed on input network
  if (statOnNw == FALSE) {
    if (!is.null(status.vector)) {
      status <- status.vector
    } else {
      status <- rep("s", num.all)
      status[sample(which(group == 1), size = i.num.all)] <- "i"
      if (groups == 2) {
        status[sample(which(group == 2), size = i.num.g2)] <- "i"
      }
      if (type == "SIR"  && !is.null(type)) {
        status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
        if (groups == 2) {
          status[sample(which(group == 2 & status == "s"),
                        size = r.num.g2)] <- "r"
        }
      }
    }
    dat <- set_attr(dat, "status", status)
  } else {
    status <- get_vertex_attribute(dat$nw[[1]], "status")
    dat <- set_attr(dat, "status", status)
  }

  origin <- get_vertex_attribute(dat$nw[[1]], "origin")
  dat <- set_attr(dat, "origin", origin)


  # Set up TEA status
  if (tergmLite == FALSE) {
    if (statOnNw == FALSE) {
      dat$nw[[1]] <- set_vertex_attribute(dat$nw[[1]], "status", status)
    }
    if (isTERGM == TRUE) {
      dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                               prefix = "testatus",
                                               value = status,
                                               onset = 1,
                                               terminus = Inf)
      dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                               prefix = "global_track",
                                               value = origin,
                                               onset = 1,
                                               terminus = Inf)
    } else {
      dat$temp$nw_list[[1]] <- set_vertex_attribute(dat$temp$nw_list[[1]],
                                                    "status", status)
    }
  }


  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  if (type == "None") {
    infTime <- rep(NA, num.all)
    idsInf <- which(status == "i")
    infTime[idsInf] <- 1
    #dat <- set_attr(dat, "infTime", infTime)
    } else {
      idsInf <- which(status == "i")
      infTime <- rep(NA, length(status))
      infTime.vector <- get_init(dat, "infTime.vector",
                                 override.null.error = TRUE)
      if (!is.null(infTime.vector)){
        infTime <- infTime.vector
        } else {
          infTime[idsInf] <- ssample(1:(-nsteps + 2),
                                     length(idsInf), replace = TRUE)
        }
      }
    dat <- set_attr(dat, "infTime", infTime)


    #min.hiv.time <- round(6.4 + 6.4)
    dat <- set_attr(dat, "stage", rep(NA, num.all))
    dat <- set_attr(dat, "stage.time", rep(NA, num.all))
    dat <- set_attr(dat, "aids.time", rep(NA, num.all))
    stage <- get_attr(dat, "stage")
    stage[idsInf] <- 0
    dat <- set_attr(dat, "stage", stage)
    dat <- track_stages(dat, at = 0, 0, idsInf)
    stage.time <- get_attr(dat, "stage.time")
    #stage.time[idsInf] <- infTime[idsInf] - min.hiv.time
    stage.time[idsInf] <- infTime[idsInf]
    dat <- set_attr(dat, "stage.time", stage.time)


    dat <- set_attr(dat, "diag.stage", rep(NA, num.all))
    diag.stage <- get_attr(dat, "diag.stage")
    diag.stage[idsInf] <- stage[idsInf]
    dat <- set_attr(dat, "diag.stage", diag.stage)

    #hiv.test.rate <- get_param(dat, "hiv.test.rate")
    hiv.test.rate.df <- get_param(dat, "hiv.test.rate")
    hiv.test.rate <- mean(hiv.test.rate.df$perc_per_day)
    dat <- set_attr(dat, "diag.time", rep(NA, num.all))
    diag.time <- get_attr(dat, "diag.time")
    #diag.time[idsInf] <- infTime[idsInf] + round(1/dat$param$hiv.test.rate)
    diag.time[idsInf] <- infTime[idsInf] + round(1/hiv.test.rate)
    dat <- set_attr(dat, "diag.time", diag.time)

    dat <- set_attr(dat, "last.neg.test",rep(NA, num.all))

    dat <- set_attr(dat, "tx.status",rep(NA, num.all))
    tx.status <- get_attr(dat, "tx.status")
    tx.status[idsInf] <- 0
    dat <- set_attr(dat, "tx.status", tx.status)

    dat <- set_attr(dat, "cuml.time.on.tx", rep(0, num.all))
    cuml.time.on.tx <- get_attr(dat, "cuml.time.on.tx")
    cuml.time.on.tx[idsInf] <- 0
    dat <- set_attr(dat, "cuml.time.on.tx", cuml.time.on.tx)


    dat <- set_attr(dat, "cuml.time.off.tx", rep(NA, num.all))
    cuml.time.off.tx <- get_attr(dat, "cuml.time.off.tx")
    cuml.time.off.tx[idsInf] <- 0
    dat <- set_attr(dat, "cuml.time.off.tx", cuml.time.off.tx)

    dat <- set_attr(dat, "tx.period.first", rep(NA, num.all))
    dat <- set_attr(dat, "tx.period.last", rep(NA, num.all))
    dat <- set_attr(dat, "tx.init.time", rep(NA, num.all))

    dat <- set_attr(dat, "count.trans", rep(0, num.all))
    #dat <- set_attr(dat, "num.neg.tests", rep(0, num.all))


  return(dat)
}
