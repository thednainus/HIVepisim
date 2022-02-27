#' @title Depature Module
#'
#' @description Module function for simulating both general and disease-related
#'              departures, including deaths, among population members.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and AIDS-mortality rate.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#'
#' @keywords module msm
#' @export
#'
departure_mig <- function(dat, at) {

  ## General departures
  #browser()
  init_sim_date <- get_param(dat, "init_date")
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  origin <- get_attr(dat, "origin")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")
  asmr <- get_param(dat, "asmr")
  rates.aids <- get_param(dat, "aids.mr")

  idsDpt.pop1 <- NULL
  idsDpt.aids.pop1 <- NULL

  idsDpt.pop2 <- NULL
  idsDpt.aids.pop2 <- NULL

  #browser()


  if(class(asmr) != "numeric"){
    # get general mortality rate by time step
    #get a2.rate for time step "at"
    dep_times <- asmr$dr_times
    #linear interpolation for time step at
    rates.all <- apply(asmr$dep_vec, 1, get_rate, init_date = init_sim_date, times = dep_times, at = at)
  } else{
    rates.all <- asmr
  }


  # Departures (not HIV-related: population 1) --------------------------------------------------
  nDepartures.pop1 <- 0
  idsElig.pop1 <- which(active == 1 & origin == "region")
  nElig.pop1 <- length(idsElig.pop1)

  whole_ages_of_elig1 <- pmin(ceiling(age[idsElig.pop1]), 80)
  drates_of_elig1 <- rates.all[whole_ages_of_elig1 - 18 + 1] #to use the correct index

  if (nElig.pop1 > 0) {
    vecDepartures.pop1 <- which(rbinom(nElig.pop1, 1, drates_of_elig1) == 1)
    if (length(vecDepartures.pop1) > 0) {
      #browser()
      idsDpt.pop1 <- idsElig.pop1[vecDepartures.pop1]
      if(any(origin[idsDpt.pop1] != "region")){
        stop("Something is wrong as an individual from global died instead of region")
      }
      nDepartures.pop1 <- length(idsDpt.pop1)
      if("i" %in% dat$attr$status[idsDpt.pop1]){
        #browser()
      }
      active[idsDpt.pop1] <- 0
      exitTime[idsDpt.pop1] <- at
    }
  }

  #browser()

  # Departures (not HIV-related: population 2) --------------------------------------------------
  nDepartures.pop2 <- 0
  idsElig.pop2 <- which(active == 1 & origin == "global")
  nElig.pop2 <- length(idsElig.pop2)

  whole_ages_of_elig2 <- pmin(ceiling(age[idsElig.pop2]), 80)
  drates_of_elig2 <- rates.all[whole_ages_of_elig2 - 18 + 1] #to use the correct index

  if (nElig.pop2 > 0) {
    vecDepartures.pop2 <- which(rbinom(nElig.pop2, 1, drates_of_elig2) == 1)
    if (length(vecDepartures.pop2) > 0) {
      #browser()
      idsDpt.pop2 <- idsElig.pop2[vecDepartures.pop2]
      if(any(origin[idsDpt.pop2] != "global")){
        stop("Something is wrong as an individual from region died instead of global")
      }
      nDepartures.pop2 <- length(idsDpt.pop2)
      active[idsDpt.pop2] <- 0
      exitTime[idsDpt.pop2] <- at
      if("i" %in% dat$attr$status[idsDpt.pop2]){
        #browser()
      }
    }
  }


  # AIDS-related departures (pop1) -----------------------------------------------------
  nDepartures.aids.pop1 <- 0
  idsElig.aids.pop1 <- which(active == 1 & stage == 4 & origin == "region")
  nElig.aids.pop1 <- length(idsElig.aids.pop1)
  if (nElig.aids.pop1 > 0) {
    vecDepartures.aids.pop1 <- which(rbinom(nElig.aids.pop1, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop1) > 0) {
      #browser()
      idsDpt.aids.pop1 <- idsElig.aids.pop1[vecDepartures.aids.pop1]
      if(any(origin[idsDpt.aids.pop1] != "region")){
        stop("Something is wrong as an individual from global died instead of region")
      }
      nDepartures.aids.pop1 <- length(idsDpt.aids.pop1)
      active[idsDpt.aids.pop1] <- 0
      exitTime[idsDpt.aids.pop1] <- at
    }
  }

  # AIDS-related departures (pop2) -----------------------------------------------------
  nDepartures.aids.pop2 <- 0
  idsElig.aids.pop2 <- which(active == 1 & stage == 4 & origin == "global")
  nElig.aids.pop2 <- length(idsElig.aids.pop2)
  if (nElig.aids.pop2 > 0) {
    vecDepartures.aids.pop2 <- which(rbinom(nElig.aids.pop2, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop2) > 0) {
      #browser()
      idsDpt.aids.pop2 <- idsElig.aids.pop2[vecDepartures.aids.pop2]
      if(any(origin[idsDpt.aids.pop2] != "global")){
        stop("Something is wrong as an individual from region died instead of global")
      }
      nDepartures.aids.pop2 <- length(idsDpt.aids.pop2)
      active[idsDpt.aids.pop2] <- 0
      exitTime[idsDpt.aids.pop2] <- at
    }
  }


  # counting the number of departures by natural causes
  # that was HIV positive
  idsDepAll.pop1 <- unique(c(idsDpt.pop1, idsDpt.aids.pop1))
  depHIV.pop1 <- intersect(idsDepAll.pop1, which(status == "i"))
  ndepHIV.pop1 <- length(depHIV.pop1)


  idsDepAll.pop2 <- unique(c(idsDpt.pop2, idsDpt.aids.pop2))
  depHIV.pop2 <- intersect(idsDepAll.pop2, which(status == "i"))
  ndepHIV.pop2 <- length(depHIV.pop2)


  #Cumulative R0 calculations
  if (length(depHIV.pop1) > 0) {
    #browser()

    newR0_pop1 <- dat$attr$count.trans[depHIV.pop1]
    r0pop1 <- data.frame(time = at, r0 = newR0_pop1)

    if(!is.null(dat$stats$R0_pop1) == TRUE){
      r0pop1 <- rbind(dat$stats$R0_pop1, r0pop1)
    }
    dat$stats$R0_pop1 <- r0pop1

  }
  if (length(depHIV.pop2) > 0) {

    #browser()

    newR0_pop2 <- dat$attr$count.trans[depHIV.pop2]
    r0pop2 <- data.frame(time = at, r0 = newR0_pop2)

    if(!is.null(dat$stats$R0_pop2) == TRUE){
      r0pop2 <- rbind(dat$stats$R0_pop2, r0pop2)
    }
    dat$stats$R0_pop2 <- r0pop2
  }


  # Output ------------------------------------------------------------------

  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  dat <- set_epi(dat, "dall_pop1.flow", at, nDepartures.pop1)
  dat <- set_epi(dat, "daids_pop1.flow", at, nDepartures.aids.pop1)
  dat <- set_epi(dat, "dhiv_pop1.flow", at, ndepHIV.pop1)

  dat <- set_epi(dat, "dall_pop2.flow", at, nDepartures.pop2)
  dat <- set_epi(dat, "daids_pop2.flow", at, nDepartures.aids.pop2)
  dat <- set_epi(dat, "dhiv_pop2.flow", at, ndepHIV.pop2)

  return(dat)
}

#' @title Depature Module
#'
#' @description Module function for simulating both general and disease-related
#'              departures, including deaths, among population members.
#'              In this function we use two departure rates: one for population
#'              1 and another for population 2.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and AIDS-mortality rate.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#'
#' @keywords module msm
#' @export
#'
departure_mig2 <- function(dat, at) {

  # attributes
  init_sim_date <- get_param(dat, "init_date")
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  origin <- get_attr(dat, "origin")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")

  #parameters
  asmr_pop1 <- get_param(dat, "asmr_pop1")
  asmr_pop2 <- get_param(dat, "asmr_pop2")
  rates.aids <- get_param(dat, "aids.mr")

  idsDpt.pop1 <- NULL
  idsDpt.aids.pop1 <- NULL

  idsDpt.pop2 <- NULL
  idsDpt.aids.pop2 <- NULL

  #browser()


  if(class(asmr_pop1) != "numeric"){
    # get general mortality rate by time step

    dep_times_pop1 <- asmr_pop1$dr_times
    dep_times_pop2 <- asmr_pop2$dr_times
    #linear interpolation for time step at
    rates.all_pop1 <- apply(asmr_pop1$dep_vec, 1, get_rate,
                            init_date = init_sim_date, times = dep_times_pop1,
                            at = at)

    rates.all_pop2 <- apply(asmr_pop2$dep_vec, 1, get_rate,
                            init_date = init_sim_date, times = dep_times_pop2,
                            at = at)
  } else{
    rates.all_pop1 <- asmr_pop1
    rates.all_pop2 <- asmr_pop2
  }


  # Departures (not HIV-related: population 1) --------------------------------------------------
  nDepartures.pop1 <- 0
  idsElig.pop1 <- which(active == 1 & origin == "region")
  nElig.pop1 <- length(idsElig.pop1)

  whole_ages_of_elig1 <- pmin(ceiling(age[idsElig.pop1]), 80)
  drates_of_elig1 <- rates.all_pop1[whole_ages_of_elig1 - 18 + 1] #to use the correct index

  if (nElig.pop1 > 0) {
    vecDepartures.pop1 <- which(rbinom(nElig.pop1, 1, drates_of_elig1) == 1)
    if((at > (10 * 365)) & (length(vecDepartures.pop1) > 0)){
      #browser()
    }
    if (length(vecDepartures.pop1) > 0) {
      #browser()
      idsDpt.pop1 <- idsElig.pop1[vecDepartures.pop1]
      if(any(origin[idsDpt.pop1] != "region")){
        stop("Something is wrong as an individual from global died instead of region")
      }
      nDepartures.pop1 <- length(idsDpt.pop1)
      if("i" %in% dat$attr$status[idsDpt.pop1]){
        #browser()
      }
      active[idsDpt.pop1] <- 0
      exitTime[idsDpt.pop1] <- at
    }
  }

  #browser()

  # Departures (not HIV-related: population 2) --------------------------------------------------
  nDepartures.pop2 <- 0
  idsElig.pop2 <- which(active == 1 & origin == "global")
  nElig.pop2 <- length(idsElig.pop2)

  whole_ages_of_elig2 <- pmin(ceiling(age[idsElig.pop2]), 80)
  drates_of_elig2 <- rates.all_pop2[whole_ages_of_elig2 - 18 + 1] #to use the correct index

  if (nElig.pop2 > 0) {
    vecDepartures.pop2 <- which(rbinom(nElig.pop2, 1, drates_of_elig2) == 1)
    if((at > (10 * 365)) & (length(vecDepartures.pop2) > 0)){
      #browser()
    }
    if (length(vecDepartures.pop2) > 0) {
      #browser()
      idsDpt.pop2 <- idsElig.pop2[vecDepartures.pop2]
      if(any(origin[idsDpt.pop2] != "global")){
        stop("Something is wrong as an individual from region died instead of global")
      }
      nDepartures.pop2 <- length(idsDpt.pop2)
      active[idsDpt.pop2] <- 0
      exitTime[idsDpt.pop2] <- at
      if("i" %in% dat$attr$status[idsDpt.pop2]){
        #browser()
      }
    }
  }


  # AIDS-related departures (pop1) -----------------------------------------------------
  nDepartures.aids.pop1 <- 0
  idsElig.aids.pop1 <- which(active == 1 & stage == 4 & origin == "region")
  nElig.aids.pop1 <- length(idsElig.aids.pop1)
  if (nElig.aids.pop1 > 0) {
    vecDepartures.aids.pop1 <- which(rbinom(nElig.aids.pop1, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop1) > 0) {
      #browser()
      idsDpt.aids.pop1 <- idsElig.aids.pop1[vecDepartures.aids.pop1]
      if(any(origin[idsDpt.aids.pop1] != "region")){
        stop("Something is wrong as an individual from global died instead of region")
      }
      nDepartures.aids.pop1 <- length(idsDpt.aids.pop1)
      active[idsDpt.aids.pop1] <- 0
      exitTime[idsDpt.aids.pop1] <- at
    }
  }

  # AIDS-related departures (pop2) -----------------------------------------------------
  nDepartures.aids.pop2 <- 0
  idsElig.aids.pop2 <- which(active == 1 & stage == 4 & origin == "global")
  nElig.aids.pop2 <- length(idsElig.aids.pop2)
  if (nElig.aids.pop2 > 0) {
    vecDepartures.aids.pop2 <- which(rbinom(nElig.aids.pop2, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop2) > 0) {
      #browser()
      idsDpt.aids.pop2 <- idsElig.aids.pop2[vecDepartures.aids.pop2]
      if(any(origin[idsDpt.aids.pop2] != "global")){
        stop("Something is wrong as an individual from region died instead of global")
      }
      nDepartures.aids.pop2 <- length(idsDpt.aids.pop2)
      active[idsDpt.aids.pop2] <- 0
      exitTime[idsDpt.aids.pop2] <- at
    }
  }


  # counting the number of departures by natural causes
  # that was HIV positive
  idsDepAll.pop1 <- unique(c(idsDpt.pop1, idsDpt.aids.pop1))
  depHIV.pop1 <- intersect(idsDepAll.pop1, which(status == "i"))
  ndepHIV.pop1 <- length(depHIV.pop1)


  idsDepAll.pop2 <- unique(c(idsDpt.pop2, idsDpt.aids.pop2))
  depHIV.pop2 <- intersect(idsDepAll.pop2, which(status == "i"))
  ndepHIV.pop2 <- length(depHIV.pop2)


  #Cumulative R0 calculations
  if (length(depHIV.pop1) > 0) {
    #browser()

    newR0_pop1 <- dat$attr$count.trans[depHIV.pop1]
    r0pop1 <- data.frame(time = at, r0 = newR0_pop1)

    if(!is.null(dat$stats$R0_pop1) == TRUE){
      r0pop1 <- rbind(dat$stats$R0_pop1, r0pop1)
    }
    dat$stats$R0_pop1 <- r0pop1

  }
  if (length(depHIV.pop2) > 0) {

    #browser()

    newR0_pop2 <- dat$attr$count.trans[depHIV.pop2]
    r0pop2 <- data.frame(time = at, r0 = newR0_pop2)

    if(!is.null(dat$stats$R0_pop2) == TRUE){
      r0pop2 <- rbind(dat$stats$R0_pop2, r0pop2)
    }
    dat$stats$R0_pop2 <- r0pop2
  }


  # Output ------------------------------------------------------------------

  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  dat <- set_epi(dat, "dall_pop1.flow", at, nDepartures.pop1)
  dat <- set_epi(dat, "daids_pop1.flow", at, nDepartures.aids.pop1)
  dat <- set_epi(dat, "dhiv_pop1.flow", at, ndepHIV.pop1)

  dat <- set_epi(dat, "dall_pop2.flow", at, nDepartures.pop2)
  dat <- set_epi(dat, "daids_pop2.flow", at, nDepartures.aids.pop2)
  dat <- set_epi(dat, "dhiv_pop2.flow", at, ndepHIV.pop2)

  return(dat)
}

#' @title Depature Module
#'
#' @description Module function for simulating both general and disease-related
#'              departures, including deaths, among population members.
#'              In this function we use two departure rates: one for population
#'              1 and another for population 2.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and AIDS-mortality rate.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#'
#' @keywords module msm
#' @export
#'
departure_mig3 <- function(dat, at) {

  # attributes
  init_sim_date <- get_param(dat, "init_date")
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  origin <- get_attr(dat, "origin")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")

  #parameters
  asmr_pop1 <- get_param(dat, "asmr_pop1")
  #asmr_pop2 <- get_param(dat, "asmr_pop2")
  rates.aids <- get_param(dat, "aids.mr")

  idsDpt.pop1 <- NULL
  idsDpt.aids.pop1 <- NULL

  idsDpt.pop2 <- NULL
  idsDpt.aids.pop2 <- NULL

  #browser()

  index <- at - 1
  n_pop1 <- get_epi(dat, "num.pop1", index)
  n_pop2 <- get_epi(dat, "num.pop2", index)

  #browser()


  if(class(asmr_pop1) != "numeric"){
    # get general mortality rate by time step

    dep_times_pop1 <- asmr_pop1$dr_times
    #dep_times_pop2 <- asmr_pop2$dr_times
    #linear interpolation for time step at
    rates.all_pop1 <- apply(asmr_pop1$dep_vec, 1, get_rate,
                            init_date = init_sim_date, times = dep_times_pop1,
                            at = at)

    #rates.all_pop2 <- apply(asmr_pop2$dep_vec, 1, get_rate,
    #                        init_date = init_sim_date, times = dep_times_pop2,
    #                        at = at)

    #to keep the arrivals balanced with arrivals in population 1
    rates.all_pop2 <- (rates.all_pop1 * n_pop1)/n_pop2

  } else{
    rates.all_pop1 <- asmr_pop1
    rates.all_pop2 <- (rates.all_pop1 * n_pop1)/n_pop2
  }


  # Departures (not HIV-related: population 1) --------------------------------------------------
  nDepartures.pop1 <- 0
  idsElig.pop1 <- which(active == 1 & origin == "region")
  nElig.pop1 <- length(idsElig.pop1)

  whole_ages_of_elig1 <- pmin(ceiling(age[idsElig.pop1]), 80)
  drates_of_elig1 <- rates.all_pop1[whole_ages_of_elig1 - 18 + 1] #to use the correct index

  if (nElig.pop1 > 0) {
    vecDepartures.pop1 <- which(rbinom(nElig.pop1, 1, drates_of_elig1) == 1)
    if((at > (10 * 365)) & (length(vecDepartures.pop1) > 0)){
      #browser()
    }
    if (length(vecDepartures.pop1) > 0) {
      #browser()
      idsDpt.pop1 <- idsElig.pop1[vecDepartures.pop1]
      if(any(origin[idsDpt.pop1] != "region")){
        stop("Something is wrong as an individual from global died instead of region")
      }
      nDepartures.pop1 <- length(idsDpt.pop1)
      if("i" %in% dat$attr$status[idsDpt.pop1]){
        #browser()
      }
      active[idsDpt.pop1] <- 0
      exitTime[idsDpt.pop1] <- at
    }
  }

  #browser()

  # Departures (not HIV-related: population 2) --------------------------------------------------
  nDepartures.pop2 <- 0
  idsElig.pop2 <- which(active == 1 & origin == "global")
  nElig.pop2 <- length(idsElig.pop2)

  whole_ages_of_elig2 <- pmin(ceiling(age[idsElig.pop2]), 80)
  drates_of_elig2 <- rates.all_pop2[whole_ages_of_elig2 - 18 + 1] #to use the correct index

  if (nElig.pop2 > 0) {
    vecDepartures.pop2 <- which(rbinom(nElig.pop2, 1, drates_of_elig2) == 1)
    if((at > (10 * 365)) & (length(vecDepartures.pop2) > 0)){
      #browser()
    }
    if (length(vecDepartures.pop2) > 0) {
      #browser()
      idsDpt.pop2 <- idsElig.pop2[vecDepartures.pop2]
      if(any(origin[idsDpt.pop2] != "global")){
        stop("Something is wrong as an individual from region died instead of global")
      }
      nDepartures.pop2 <- length(idsDpt.pop2)
      active[idsDpt.pop2] <- 0
      exitTime[idsDpt.pop2] <- at
      if("i" %in% dat$attr$status[idsDpt.pop2]){
        #browser()
      }
    }
  }


  # AIDS-related departures (pop1) -----------------------------------------------------
  nDepartures.aids.pop1 <- 0
  idsElig.aids.pop1 <- which(active == 1 & stage == 4 & origin == "region")
  nElig.aids.pop1 <- length(idsElig.aids.pop1)
  if (nElig.aids.pop1 > 0) {
    vecDepartures.aids.pop1 <- which(rbinom(nElig.aids.pop1, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop1) > 0) {
      #browser()
      idsDpt.aids.pop1 <- idsElig.aids.pop1[vecDepartures.aids.pop1]
      if(any(origin[idsDpt.aids.pop1] != "region")){
        stop("Something is wrong as an individual from global died instead of region")
      }
      nDepartures.aids.pop1 <- length(idsDpt.aids.pop1)
      active[idsDpt.aids.pop1] <- 0
      exitTime[idsDpt.aids.pop1] <- at
    }
  }

  # AIDS-related departures (pop2) -----------------------------------------------------
  nDepartures.aids.pop2 <- 0
  idsElig.aids.pop2 <- which(active == 1 & stage == 4 & origin == "global")
  nElig.aids.pop2 <- length(idsElig.aids.pop2)
  if (nElig.aids.pop2 > 0) {
    vecDepartures.aids.pop2 <- which(rbinom(nElig.aids.pop2, 1, rates.aids) == 1)
    if (length(vecDepartures.aids.pop2) > 0) {
      #browser()
      idsDpt.aids.pop2 <- idsElig.aids.pop2[vecDepartures.aids.pop2]
      if(any(origin[idsDpt.aids.pop2] != "global")){
        stop("Something is wrong as an individual from region died instead of global")
      }
      nDepartures.aids.pop2 <- length(idsDpt.aids.pop2)
      active[idsDpt.aids.pop2] <- 0
      exitTime[idsDpt.aids.pop2] <- at
    }
  }


  # counting the number of departures by natural causes
  # that was HIV positive
  idsDepAll.pop1 <- unique(c(idsDpt.pop1, idsDpt.aids.pop1))
  depHIV.pop1 <- intersect(idsDepAll.pop1, which(status == "i"))
  ndepHIV.pop1 <- length(depHIV.pop1)


  idsDepAll.pop2 <- unique(c(idsDpt.pop2, idsDpt.aids.pop2))
  depHIV.pop2 <- intersect(idsDepAll.pop2, which(status == "i"))
  ndepHIV.pop2 <- length(depHIV.pop2)


  #Cumulative R0 calculations
  if (length(depHIV.pop1) > 0) {
    #browser()

    newR0_pop1 <- dat$attr$count.trans[depHIV.pop1]
    r0pop1 <- data.frame(time = at, r0 = newR0_pop1)

    if(!is.null(dat$stats$R0_pop1) == TRUE){
      r0pop1 <- rbind(dat$stats$R0_pop1, r0pop1)
    }
    dat$stats$R0_pop1 <- r0pop1

  }
  if (length(depHIV.pop2) > 0) {

    #browser()

    newR0_pop2 <- dat$attr$count.trans[depHIV.pop2]
    r0pop2 <- data.frame(time = at, r0 = newR0_pop2)

    if(!is.null(dat$stats$R0_pop2) == TRUE){
      r0pop2 <- rbind(dat$stats$R0_pop2, r0pop2)
    }
    dat$stats$R0_pop2 <- r0pop2
  }


  # Output ------------------------------------------------------------------

  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  dat <- set_epi(dat, "dall_pop1.flow", at, nDepartures.pop1)
  dat <- set_epi(dat, "daids_pop1.flow", at, nDepartures.aids.pop1)
  dat <- set_epi(dat, "dhiv_pop1.flow", at, ndepHIV.pop1)

  dat <- set_epi(dat, "dall_pop2.flow", at, nDepartures.pop2)
  dat <- set_epi(dat, "daids_pop2.flow", at, nDepartures.aids.pop2)
  dat <- set_epi(dat, "dhiv_pop2.flow", at, ndepHIV.pop2)

  return(dat)
}
