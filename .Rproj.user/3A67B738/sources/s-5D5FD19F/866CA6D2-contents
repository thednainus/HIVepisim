#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams EpiModelHIV::aging_msm
#'
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes. Summary statistics for
#' number of new diagnoses, total number of tests, and total number of negative
#' tests are calculated and stored on \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
hivtest_msm <- function(dat, at) {

  ## Variables

  # Attributes -----
  diag.status <- get_attr(dat, "diag.status")
  diag.time <- get_attr(dat, "diag.time")
  diag.stage <- get_attr(dat, "diag.stage")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  stage <- get_attr(dat, "stage")
  entrTime <- get_attr(dat, "entrTime")
  last.neg.test <- get_attr(dat, "last.neg.test")
  num.neg.tests <- get_attr(dat, "num.neg.tests")

  tsincelntst <- at - last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - entrTime[is.na(tsincelntst)]

  # Parameters ----
  hiv.test.rate <- get_param(dat, "hiv.test.rate")
  twind.int <- get_param(dat, "test.window.int")


  # General interval testing
  elig <- which((diag.status == 0 | is.na(diag.status)))


    # Testing rates
  rates <- hiv.test.rate
  idsTstGen <- elig[rbinom(length(elig), 1, rates) == 1]

  tstAll <- idsTstGen

  tstPos <- tstAll[status[tstAll] == 1 & infTime[tstAll] <= at - twind.int]
  tstNeg <- setdiff(tstAll, tstPos)

  # Outputs  -----
  last.neg.test[tstNeg] <- at
  diag.status[tstPos] <- 1
  diag.time[tstPos] <- at
  diag.stage[tstPos] <- stage[tstPos]

  # Summary stats
  if (at >= 52*65) {
    num.neg.tests[tstNeg] <- num.neg.tests[tstNeg] + 1
  }

  dat <- set_attr(dat, "last.neg.test", last.neg.test)
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_attr(dat, "diag.time", diag.time)
  dat <- set_attr(dat, "diag.stage", diag.stage)


  dat <- set_epi(dat, "tot.tests", at, length(tstAll))
  dat <- set_epi(dat, "tot.neg.tests", at, length(tstNeg))
  # number of new diagnoses by timing
  dat <- set_epi(dat, "newDx", at,  length(tstPos))

  return(dat)
}

