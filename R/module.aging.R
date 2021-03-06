
#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{age}.
#'
#' @export
#'

aging_msm <- function(dat, at) {

  ## Parameters -----------------
  time.unit <- get_param(dat,"time.unit")

  ## Attributes ------------------
  age <- get_attr(dat, "age")
  active <- get_attr(dat, "active")

  ## Updates ----------------------
  age[active == 1] <- age[active == 1] + time.unit/365

  ## Set updated age attribute on dat object
  dat <- set_attr(dat, "age", age)

  return(dat)
}
