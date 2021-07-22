#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams EpiModel::arrivals.net
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
  active <- get_attr(dat, "active")
  origin <- get_attr(dat, "origin")
  diag.status <- get_attr(dat, "diag.status")
  diag.time <- get_attr(dat, "diag.time")
  diag.stage <- get_attr(dat, "diag.stage")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  stage <- get_attr(dat, "stage")
  entrTime <- get_attr(dat, "entrTime")
  last.neg.test <- get_attr(dat, "last.neg.test")
  #num.neg.tests <- get_attr(dat, "num.neg.tests")

  tsincelntst <- at - last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - entrTime[is.na(tsincelntst)]

  # Parameters ----
  hiv.test.rate <- get_param(dat, "hiv.test.rate")
  twind.int <- get_param(dat, "test.window.int")


  # General interval testing
  elig_pop1 <- which(active == 1 & origin == "region" & (diag.status == 0 | is.na(diag.status)))
  elig_pop2 <- which(active == 1 & origin == "global" & (diag.status == 0 | is.na(diag.status)))


    # Testing rates
  #browser()
  idsTstGen_pop1 <- elig_pop1[rbinom(length(elig_pop1), 1, hiv.test.rate) == 1]
  idsTstGen_pop2 <- elig_pop2[rbinom(length(elig_pop2), 1, hiv.test.rate) == 1]


  tstPos_pop1 <- idsTstGen_pop1[status[idsTstGen_pop1] == "i" & infTime[idsTstGen_pop1] <= at - twind.int]
  tstPos_pop2 <- idsTstGen_pop2[status[idsTstGen_pop2] == "i" & infTime[idsTstGen_pop2] <= at - twind.int]

  tstNeg_pop1 <- setdiff(idsTstGen_pop1, tstPos_pop1)
  tstNeg_pop2 <- setdiff(idsTstGen_pop2, tstPos_pop2)

  #set diagnosis time
  if(dat$control$save.stats == TRUE){
    #browser()
    all_pos_tes_ids <- union(tstPos_pop1, tstPos_pop2)
    if(length(all_pos_tes_ids) > 0){
      dat <- set_diagnosis_time(dat, at, union(tstPos_pop1, tstPos_pop2))
    }
  }

  # Outputs  -----
  last.neg.test[union(tstNeg_pop1, tstNeg_pop2)] <- at
  diag.status[union(tstPos_pop1, tstPos_pop2)] <- 1
  diag.time[union(tstPos_pop1, tstPos_pop2)] <- at
  diag.stage[union(tstPos_pop1, tstPos_pop2)] <- stage[union(tstPos_pop1, tstPos_pop2)]

  # Summary stats
  #if (at >= 52*65) {
  #  num.neg.tests[tstNeg] <- num.neg.tests[tstNeg] + 1
  #}

  dat <- set_attr(dat, "last.neg.test", last.neg.test)
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_attr(dat, "diag.time", diag.time)
  dat <- set_attr(dat, "diag.stage", diag.stage)


  dat <- set_epi(dat, "tot.tests_pop1", at, length(idsTstGen_pop1))
  dat <- set_epi(dat, "tot.tests_pop2", at, length(idsTstGen_pop2))
  dat <- set_epi(dat, "tot.neg.tests_pop1", at, length(tstNeg_pop1))
  dat <- set_epi(dat, "tot.neg.tests_pop2", at, length(tstNeg_pop2))
  # number of new diagnoses by timing
  dat <- set_epi(dat, "newDx_pop1", at,  length(tstPos_pop1))
  dat <- set_epi(dat, "newDx_pop2", at,  length(tstPos_pop2))

  return(dat)
}

