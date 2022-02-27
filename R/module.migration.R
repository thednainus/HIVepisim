#' @title Migration Module
#'
#' @description Module function for migrations between two populations
#'
#' @inheritParams EpiModel::arrivals.net
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


  #browser()

  if(m12.rate == 0){
    nMig12 <- 0
  }

  if(m21.rate == 0){
    nMig21 <- 0
  }



  if(m12.rate != 0){
    ## Process
    # migrations for population 1 to population 2

    nMig12 <- 0

    # Attributes ----
    active <- get_attr(dat, "active")
    origin <- get_attr(dat, "origin")
    migrant <- get_attr(dat, "migrant")
    migrationTime <- get_attr(dat, "migrationTime")


    idsMigs12 <- which(active == 1 & origin == "region")
    nidsMigs12 <- length(idsMigs12)

    ## Update Attr
    if (nidsMigs12 > 0) {
      # because this is a migration from population 1 to population 2
      # the migrants will be labeled as origin = "global" (because it
      # is not in population 1 anymore)
      #random select the number fo migrants from the network to migrate
      vecMigrations.12 <- which(rbinom(nidsMigs12, 1, m12.rate) == 1)
      if (length(vecMigrations.12) > 0) {
        idsMigs12.all <- idsMigs12[vecMigrations.12]
        #browser()
        if(any(origin[idsMigs12.all] != "region")){
          stop("Something is wrong as someone in global migrated to region, instead of
               region to global")
        }
        nMig12 <- length(idsMigs12.all)
        origin[idsMigs12.all] <- "global"
        migrant[idsMigs12.all] <- 12
        migrationTime[idsMigs12.all] <- at
        dat <- set_attr(dat, "origin", origin)
        dat <- set_attr(dat, "migrant", migrant)
        dat <- set_attr(dat, "migrationTime", migrationTime)
        dat <- track_origin(dat, at, "12", idsMigs12.all)
        #dat$nw[[1]] <- deactivate.edges(dat$nw[[1]], onset = at,
        #                                terminus = Inf,
        #                                e = get.edgeIDs(x = dat$nw[[1]], v = idsMigs12.all))
        #dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]], prefix = "global_track",
        #                                         value = 1, onset = at,
        #                                         terminus = Inf, v = idsMigs12.all)
      }
    }
  }

  if(m21.rate != 0){
    ## Process
    # migrations for population 1 to population 2
    #size of population 1

    nMig21 <- 0

    # Attributes ----
    active <- get_attr(dat, "active")
    origin <- get_attr(dat, "origin")
    migrant <- get_attr(dat, "migrant")
    migrationTime <- get_attr(dat, "migrationTime")


    idsMigs21 <- which(active == 1 & origin == "global")
    nidsMigs21 <- length(idsMigs21)

    ## Update Attr
    if (nidsMigs21 > 0) {
      # because this is a migration from population 1 to population 2
      # the migrants will be labeled as origin = "global" (because it
      # is not in population 1 anymore)
      #random select the number of migrants from the network to migrate
      vecMigrations.21 <- which(rbinom(nidsMigs21, 1, m21.rate) == 1)
      if (length(vecMigrations.21) > 0) {
        idsMigs21.all <- idsMigs21[vecMigrations.21]
        #browser()
        if(any(origin[idsMigs21.all] != "global")){
          stop("Something is wrong as someone in region migrated to global, instead of
               global to region")
        }
        nMig21 <- length(idsMigs21.all)
        origin[idsMigs21.all] <- "region"
        migrant[idsMigs21.all] <- 21
        migrationTime[idsMigs21.all] <- at
        dat <- set_attr(dat, "origin", origin)
        dat <- set_attr(dat, "migrant", migrant)
        dat <- set_attr(dat, "migrationTime", migrationTime)
        dat <- track_origin(dat, at, "21", idsMigs21.all)
        #nMig21 <- nMigrations21.all
        #dat$nw[[1]] <- deactivate.edges(dat$nw[[1]], onset = at,
        #                                terminus = Inf,
        #                                e = get.edgeIDs(x = dat$nw[[1]], v = idsMigs21.all))
      }
    }
  }


  ## Output ----
  dat <- set_epi(dat, "nArrivals_mig1", at, nMig21)
  dat <- set_epi(dat, "nArrivals_mig2", at, nMig12)


  return(dat)
}


#' @title Migration Module
#'
#' @description Module function for migrations between two populations
#'
#' @inheritParams EpiModel::arrivals.net
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
migration2 <- function(dat, at){
  #browser()
  ## Parameters ----
  # migration rate from population 1 to population 2
  m12.migrants <- get_param(dat, "m12.migrants")
  # migration rate from population 2 to population 1
  m21.migrants <- get_param(dat, "m21.migrants")


  #browser()

  if(m12.migrants == 0){
    nMig12 <- 0
  }

  if(m21.migrants == 0){
    nMig21 <- 0
  }

  index <- at - 1
  n_pop1 <- get_epi(dat, "num.pop1", index)
  n_pop2 <- get_epi(dat, "num.pop2", index)

  m12.rate = 1/((n_pop1/m12.migrants)*365)
  m21.rate = 1/((n_pop2/m21.migrants)*365)



  if(m12.rate != 0){
    ## Process
    # migrations for population 1 to population 2

    nMig12 <- 0

    # Attributes ----
    active <- get_attr(dat, "active")
    origin <- get_attr(dat, "origin")
    migrant <- get_attr(dat, "migrant")
    migrationTime <- get_attr(dat, "migrationTime")


    idsMigs12 <- which(active == 1 & origin == "region")
    nidsMigs12 <- length(idsMigs12)

    ## Update Attr
    if (nidsMigs12 > 0) {
      # because this is a migration from population 1 to population 2
      # the migrants will be labeled as origin = "global" (because it
      # is not in population 1 anymore)
      #random select the number fo migrants from the network to migrate
      vecMigrations.12 <- which(rbinom(nidsMigs12, 1, m12.rate) == 1)
      if (length(vecMigrations.12) > 0) {
        idsMigs12.all <- idsMigs12[vecMigrations.12]
        #browser()
        if(any(origin[idsMigs12.all] != "region")){
          stop("Something is wrong as someone in global migrated to region, instead of
               region to global")
        }
        nMig12 <- length(idsMigs12.all)
        origin[idsMigs12.all] <- "global"
        migrant[idsMigs12.all] <- 12
        migrationTime[idsMigs12.all] <- at
        dat <- set_attr(dat, "origin", origin)
        dat <- set_attr(dat, "migrant", migrant)
        dat <- set_attr(dat, "migrationTime", migrationTime)
        dat <- track_origin(dat, at, "12", idsMigs12.all)
        #dat$nw[[1]] <- deactivate.edges(dat$nw[[1]], onset = at,
        #                                terminus = Inf,
        #                                e = get.edgeIDs(x = dat$nw[[1]], v = idsMigs12.all))
        #dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]], prefix = "global_track",
        #                                         value = 1, onset = at,
        #                                         terminus = Inf, v = idsMigs12.all)
      }
    }
  }

  if(m21.rate != 0){
    ## Process
    # migrations for population 1 to population 2
    #size of population 1

    nMig21 <- 0

    # Attributes ----
    active <- get_attr(dat, "active")
    origin <- get_attr(dat, "origin")
    migrant <- get_attr(dat, "migrant")
    migrationTime <- get_attr(dat, "migrationTime")


    idsMigs21 <- which(active == 1 & origin == "global")
    nidsMigs21 <- length(idsMigs21)

    ## Update Attr
    if (nidsMigs21 > 0) {
      # because this is a migration from population 1 to population 2
      # the migrants will be labeled as origin = "global" (because it
      # is not in population 1 anymore)
      #random select the number of migrants from the network to migrate
      vecMigrations.21 <- which(rbinom(nidsMigs21, 1, m21.rate) == 1)
      if (length(vecMigrations.21) > 0) {
        idsMigs21.all <- idsMigs21[vecMigrations.21]
        #browser()
        if(any(origin[idsMigs21.all] != "global")){
          stop("Something is wrong as someone in region migrated to global, instead of
               global to region")
        }
        nMig21 <- length(idsMigs21.all)
        origin[idsMigs21.all] <- "region"
        migrant[idsMigs21.all] <- 21
        migrationTime[idsMigs21.all] <- at
        dat <- set_attr(dat, "origin", origin)
        dat <- set_attr(dat, "migrant", migrant)
        dat <- set_attr(dat, "migrationTime", migrationTime)
        dat <- track_origin(dat, at, "21", idsMigs21.all)
        #nMig21 <- nMigrations21.all
        #dat$nw[[1]] <- deactivate.edges(dat$nw[[1]], onset = at,
        #                                terminus = Inf,
        #                                e = get.edgeIDs(x = dat$nw[[1]], v = idsMigs21.all))
      }
    }
  }


  ## Output ----
  dat <- set_epi(dat, "nArrivals_mig1", at, nMig21)
  dat <- set_epi(dat, "nArrivals_mig2", at, nMig12)


  return(dat)
}


