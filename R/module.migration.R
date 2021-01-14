#' @title Migration Module
#'
#' @description Module function for migrations between two populations
#'
#' @inheritParams EpiModelHIV::aging_msm
#'
#' @details
#' New migrations from both directions (from population 1 to 2; and
#' from population 2 to 1) will happen based on migration rates and population
#' size for each population stocastically determined with with draws from
#' Poisson distributions. For each new entry, a set of attributes is added for
#' that node, and the nodes are added onto the network objects.
#' Only attributes that are a part of the network model
#' formulae are updated as vertex attributes on the network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm
#' @export
#'
migration <- function(dat, at){

  ## Parameters ----
  # migration rate from population 1 to population 2
  m12.rate <- get_param(dat, "m12.rate")
  # migration rate from population 2 to population 1
  m21.rate <- get_param(dat, "m21.rate")

  # Attributes ----
  origin <- get_attr(dat, "origin")
  migrant <- get_attr(dat, "migrant")
  migrationTime <- get_attr(dat, "migrationTime")


  if(m12.rate != 0){
    ## Process
    # migrations for population 1 to population 2
    #size of population 1
    pop1 <- dat$epi$num1[at - 1]

    nMigExp12 <- pop1 * m12.rate
    nMig12 <- rpois(1, nMigExp12)

    ## Update Attr
    if (nMig12 > 0) {
      # because this is a migration from population 1 to population 2
      # the migrants will be labeled as origin = "global" (because it
      # is not in population 1 anymore)
      #random select the number fo migrants from the network to migrate
      idsNewMig12 <- sample(which(origin == "region"), nMig12)
      origin[idsNewMig12] <- "global"
      migrant[idsNewMig12] <- 12
      migrationTime[idsNewMig12] <- at
      dat <- set_attr(dat, "origin", origin)
      dat <- set_attr(dat, "migrant", migrant)
      dat <- set_attr(dat, "migrationTime", migrationTime)
    }
  } else {
    nMig12 <-  0
  }

  if(m21.rate != 0){
    # migrations for population 1 to population 2
    #size of population 1
    pop2 <- dat$epi$num2[at - 1]

    # get the number of migrants
    nMigExp21 <- pop2 * m21.rate
    nMig21 <- rpois(1, nMigExp21)


    if (nMig21 > 0) {
      # because this is a migration from population 2 to population 1
      # the migrants will be labeled as origin = "region" (because it
      # is not in population 2 anymore)
      #random select the number fo migrants from the network to migrate
      idsNewMig21 <- sample(which(origin == "global"), nMig21)
      # reset the attribute origin
      origin[idsNewMig21] <- "region"
      migrant[idsNewMig21] <- 21
      migrationTime[idsNewMig21] <- at
      dat <- set_attr(dat, "origin", origin)
      dat <- set_attr(dat, "migrant", migrant)
      dat <- set_attr(dat, "migrationTime", migrationTime)
    }
  } else {
    nMig21 <- 0
  }


  ## Output ----
  dat <- set_epi(dat, "nArrivals_mig1", at, nMig21)
  dat <- set_epi(dat, "nArrivals_mig2", at, nMig12)


  return(dat)
}
