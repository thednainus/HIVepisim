
#' @title Disease Progression Module
#'
#' @description Module function for HIV disease progression
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' HIV disease is divided into four stages: stage 0 or acute and early HIV infection,
#' acute rising, and stage 1 to stage 4 (stage 4 is AIDS) as explained in
#' Cori et al. 2015.
#'
#' @return
#' This function returns the \code{dat} object after updating the disease stage
#' of infected individuals.
#'
#'
#' @references
#' Cori A, Pickles M, van Sighem A, et al. CD4+ cell dynamics in untreated
#' HIV-1 infection: overall rates, and effects of age, viral load, sex and
#' calendar time. AIDS. 2015;29(18):2435-2446.
#' doi:10.1097/QAD.0000000000000854
#'
#' @export
#'
hivprogress_msm <- function(dat, at) {

  # Attributes ----------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  cuml.time.on.tx <- get_attr(dat, "cuml.time.on.tx")
  cuml.time.off.tx <- get_attr(dat, "cuml.time.off.tx")
  stage <- get_attr(dat, "stage")
  stage.time <- get_attr(dat, "stage.time")
  aids.time <- get_attr(dat, "aids.time")
  tx.status <- get_attr(dat, "tx.status")

  # Parameters -------------
  stage_prog_rate0 <- get_param(dat, "stage_prog_rate0")
  stage_prog_rate1 <- get_param(dat, "stage_prog_rate1")
  stage_prog_rate2 <- get_param(dat, "stage_prog_rate2")
  stage_prog_rate3 <- get_param(dat, "stage_prog_rate3")


  f1 <- get_param(dat, "f1")
  f2 <- get_param(dat, "f2")
  f3 <- get_param(dat, "f3")
  f4 <- get_param(dat, "f4")


  ## Process ----------

  # Increment day
  stage.time[active == 1] <- stage.time[active == 1] + 1

  # Change from stage 0 to stage 1 of HIV infection
  # Note that this is following Cori et al. 2015
  # Stage 0 the acute and early HIV infection

  idsStage0 <- which(active == 1 & stage == 0 & tx.status == 0)
  # ids that will progress to another stage following probabilities f1,f2,f3 and f4
  newidsStage0 <- idsStage0[rbinom(length(idsStage0), 1, stage_prog_rate0) == 1]

  # new ids that will progress from stage 0 to stages 1 or stage 2 or stage 3 or stage 4
  # this is following Cori et al 2015
  stages <- sample(x = 1:4, size = length(newidsStage0),
                   replace = TRUE, prob = c(f1, f2, f3, f4))


  # Change from stage 1 to stage 2 of HIV infection following Cori et al. 2015
  idsStage1 <- which(active == 1 & stage == 1 & tx.status == 0)
  newidsStage2 <- idsStage1[rbinom(length(idsStage1), 1, stage_prog_rate1) == 1]



  # Change from stage 2 to stage 3 of HIV infection following Cori et al. 2015
  idsStage2 <- which(active == 1 & stage == 2 & tx.status == 0)
  newidsStage3 <- idsStage2[rbinom(length(idsStage2), 1, stage_prog_rate2) == 1]



  # Change stage to AIDS
  # Change from stage 3 to stage 4 (in Cori et al 2015)
  # an individual will only change to AIDS if not on treatment
  aids.tx.naive <- which(active == 1 & status == "i" & cuml.time.on.tx == 0 & stage == 3)
  aids.off.tx <- which(active == 1 & status == "i" & tx.status == 0 & cuml.time.off.tx > 0 & stage == 3)
  isAIDS <- c(aids.tx.naive, aids.off.tx)

  newidsAIDS <- isAIDS[rbinom(length(isAIDS), 1, stage_prog_rate3) == 1]

  #stage[isAIDS] <- 4
  #stage.time[isAIDS] <- 1
  #aids.time[isAIDS] <- at


  # change stage accordingly
  stage[newidsStage0] <- stages
  stage.time[newidsStage0] <- 1

  stage[newidsStage2] <- 2
  stage.time[newidsStage2] <- 1

  stage[newidsStage3] <- 3
  stage.time[newidsStage3] <- 1

  stage[newidsAIDS] <- 4
  stage.time[newidsAIDS] <- 1
  aids.time[newidsAIDS] <- at



  ## Output ------
  dat <- set_attr(dat, "stage", stage)
  dat <- set_attr(dat, "stage.time", stage.time)
  dat <- set_attr(dat, "aids.time", aids.time)

  dat <- set_epi(dat, "new.aids.tot", at, length(isAIDS))


  return(dat)
}

