
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
  origin <- get_attr(dat, "origin")
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

  # for "region" population ----


  # Change from stage 0 to stage 1 of HIV infection
  # Note that this is following Cori et al. 2015
  # Stage 0 the acute and early HIV infection

  idsStage0_pop1 <- which(active == 1 &
                            stage == 0 &
                            tx.status == 0 &
                            status == "i" &
                            origin == "region")
  # ids that will progress to another stage following probabilities f1,f2,f3 and f4
  newidsStage0_pop1 <- idsStage0_pop1[rbinom(length(idsStage0_pop1), 1,
                                             stage_prog_rate0) == 1]

  # new ids that will progress from stage 0 to stages 1 or stage 2 or stage 3 or stage 4
  # this is following Cori et al 2015
  stages_pop1 <- sample(x = 1:4, size = length(newidsStage0_pop1),
                   replace = TRUE, prob = c(f1, f2, f3, f4))


  # Change from stage 1 to stage 2 of HIV infection following Cori et al. 2015
  idsStage1_pop1 <- which(active == 1 &
                            stage == 1 &
                            tx.status == 0 &
                            status == "i" &
                            origin == "region")
  newidsStage2_pop1 <- idsStage1_pop1[rbinom(length(idsStage1_pop1), 1,
                                             stage_prog_rate1) == 1]



  # Change from stage 2 to stage 3 of HIV infection following Cori et al. 2015
  idsStage2_pop1 <- which(active == 1 &
                            stage == 2 &
                            tx.status == 0 &
                            status == "i" &
                            origin == "region")
  newidsStage3_pop1 <- idsStage2_pop1[rbinom(length(idsStage2_pop1), 1,
                                             stage_prog_rate2) == 1]



  # Change stage to AIDS
  # Change from stage 3 to stage 4 (in Cori et al 2015)
  # an individual will only change to AIDS if not on treatment
  aids.tx.naive_pop1 <- which(active == 1 &
                                status == "i" &
                                cuml.time.on.tx == 0 &
                                stage == 3 &
                                origin == "region")
  aids.off.tx_pop1 <- which(active == 1 & status == "i" &
                              tx.status == 0 &
                              cuml.time.off.tx > 0 &
                              stage == 3 &
                              origin == "region")
  isAIDS_pop1 <- c(aids.tx.naive_pop1, aids.off.tx_pop1)

  newidsAIDS_pop1 <- isAIDS_pop1[rbinom(length(isAIDS_pop1), 1,
                                        stage_prog_rate3) == 1]


  # for "global" population ----

  # Change from stage 0 to stage 1 of HIV infection
  # Note that this is following Cori et al. 2015
  # Stage 0 the acute and early HIV infection

  idsStage0_pop2 <- which(active == 1 &
                            stage == 0 &
                            tx.status == 0 &
                            status == "i" &
                            origin == "global")
  # ids that will progress to another stage following probabilities f1,f2,f3 and f4
  newidsStage0_pop2 <- idsStage0_pop2[rbinom(length(idsStage0_pop2), 1,
                                             stage_prog_rate0) == 1]

  # new ids that will progress from stage 0 to stages 1 or stage 2 or stage 3 or stage 4
  # this is following Cori et al 2015
  stages_pop2 <- sample(x = 1:4, size = length(newidsStage0_pop2),
                        replace = TRUE, prob = c(f1, f2, f3, f4))


  # Change from stage 1 to stage 2 of HIV infection following Cori et al. 2015
  idsStage1_pop2 <- which(active == 1 &
                            stage == 1 &
                            tx.status == 0 &
                            status == "i" &
                            origin == "global")
  newidsStage2_pop2 <- idsStage1_pop2[rbinom(length(idsStage1_pop2), 1,
                                             stage_prog_rate1) == 1]



  # Change from stage 2 to stage 3 of HIV infection following Cori et al. 2015
  idsStage2_pop2 <- which(active == 1 &
                            stage == 2 &
                            tx.status == 0 &
                            status == "i" &
                            origin == "global")
  newidsStage3_pop2 <- idsStage2_pop2[rbinom(length(idsStage2_pop2), 1,
                                             stage_prog_rate2) == 1]



  # Change stage to AIDS
  # Change from stage 3 to stage 4 (in Cori et al 2015)
  # an individual will only change to AIDS if not on treatment
  aids.tx.naive_pop2 <- which(active == 1 &
                                status == "i" &
                                cuml.time.on.tx == 0 &
                                stage == 3 &
                                origin == "global")
  aids.off.tx_pop2 <- which(active == 1 &
                              status == "i" &
                              tx.status == 0 &
                              cuml.time.off.tx > 0 &
                              stage == 3 &
                              origin == "global")
  isAIDS_pop2 <- c(aids.tx.naive_pop2, aids.off.tx_pop2)

  newidsAIDS_pop2 <- isAIDS_pop2[rbinom(length(isAIDS_pop2), 1,
                                        stage_prog_rate3) == 1]


  # change stage accordingly : "region" population ----
  if(length(newidsStage0_pop1) > 0){
    stage[newidsStage0_pop1] <- stages_pop1
    stage.time[newidsStage0_pop1] <- 1
  }

  if(length(newidsStage2_pop1) > 0){
    stage[newidsStage2_pop1] <- 2
    stage.time[newidsStage2_pop1] <- 1
  }

  if(length(newidsStage3_pop1) > 0){
    stage[newidsStage3_pop1] <- 3
    stage.time[newidsStage3_pop1] <- 1
  }

  if(length(newidsAIDS_pop1) > 0){
    stage[newidsAIDS_pop1] <- 4
    stage.time[newidsAIDS_pop1] <- 1
    aids.time[newidsAIDS_pop1] <- at
  }


  # change stage accordingly : "global" population ----
  if(length(newidsStage0_pop2) > 0){
    stage[newidsStage0_pop2] <- stages_pop2
    stage.time[newidsStage0_pop2] <- 1
  }

  if(length(newidsStage2_pop2) > 0){
    stage[newidsStage2_pop2] <- 2
    stage.time[newidsStage2_pop2] <- 1
  }

  if(length(newidsStage3_pop2) > 0){
    stage[newidsStage3_pop2] <- 3
    stage.time[newidsStage3_pop2] <- 1
  }

  if(length(newidsAIDS_pop2) > 0){
    stage[newidsAIDS_pop2] <- 4
    stage.time[newidsAIDS_pop2] <- 1
    aids.time[newidsAIDS_pop2] <- at
  }



  #browser()
    #track HIV stage of infection
  if(dat$control$save.stats == TRUE){
    #browser()

    if(length(newidsStage0_pop1) > 0){
      #browser()
      dat <- track_stages(dat, at, stages_pop1, newidsStage0_pop1)
    }

    if(length(newidsStage0_pop2) > 0){
      #browser()
      dat <- track_stages(dat, at, stages_pop2, newidsStage0_pop2)
    }


    if(length(newidsStage2_pop1) > 0 | length(newidsStage2_pop2) > 0){
      #browser()
      dat <- track_stages(dat, at, 2, c(newidsStage2_pop1, newidsStage2_pop2))
    }

    if(length(newidsStage3_pop1) > 0 | length(newidsStage3_pop2) > 0){
      #browser()
      dat <- track_stages(dat, at, 3, c(newidsStage3_pop1, newidsStage3_pop2))
    }

    if(length(newidsAIDS_pop1) > 0 | length(newidsAIDS_pop2) > 0){
      #browser()
      dat <- track_stages(dat, at, 4, c(newidsAIDS_pop1, newidsAIDS_pop2))
    }

  }

  ## Output ------
  dat <- set_attr(dat, "stage", stage)
  dat <- set_attr(dat, "stage.time", stage.time)
  dat <- set_attr(dat, "aids.time", aids.time)

  #dat <- set_epi(dat, "new.aids.tot", at, length(isAIDS))


  return(dat)
}

