#' Function to append core attributes
#'
#'
#' @inheritParams EpiModel::append_core_attr
#'
#' @return
#'
#' @details
#' This function is based on \link[EpiModel]{append_core_attr}. I added migrationTime
#' as another core attribute.
#' @export
#'
#'
append_core_attr_mig <-  function (dat, at, n.new)
{
  dat <- append_attr(dat, "active", 1, n.new)
  dat <- append_attr(dat, "entrTime", at, n.new)
  dat <- append_attr(dat, "exitTime", NA, n.new)
  dat <- append_attr(dat, "migrationTime", NA, n.new)
  #dat <- update_uids(dat, n.new)
  dat <- update_unique_ids(dat, n.new)
  return(dat)
}


#' @title Create the uids for the new nodes
#'
#' @description This function is called by `append_core_attr` and append new
#' uids to the created nodes. It also keeps track of the already used uids with
#' the /code{dat[["_last_uid"]]} variable
#'
#' @param dat a Master list object of network models
#' @param n.new the number of new nodes to give \code{uid} to
#'
#' @return the Master list object of network models (\code{dat})
#'
#' @keywords internal
# update_uids <- function(dat, n.new) {
#   last_uid <- if (is.null(dat[["_last_uid"]])) 0L else dat[["_last_uid"]]
#   next_uids <- seq_len(n.new) + last_uid
#   dat[["_last_uid"]] <- last_uid + as.integer(n.new)
#   dat <- append_attr(dat, "uid", next_uids, n.new)
#
#   return(dat)
# }


#' @title Create the unique_ids for the new nodes
#'
#' @description This function is called by `append_core_attr` and append new
#' unique_ids to the created nodes. It also keeps track of the already used
#' unique_ids with the /code{dat[["_last_unique_id"]]} variable
#'
#' @param dat a Master list object of network models
#' @param n.new the number of new nodes to give \code{unique_id} to
#'
#' @return the Master list object of network models (\code{dat})
#'
#' @keywords internal
#' @export
update_unique_ids <- function(dat, n.new) {
  last_unique_id <- if (is.null(dat[["_last_unique_id"]])) 0L
  else dat[["_last_unique_id"]]
  next_unique_ids <- seq_len(n.new) + last_unique_id
  dat[["_last_unique_id"]] <- last_unique_id + as.integer(n.new)
  dat <- append_attr(dat, "unique_id", next_unique_ids, n.new)

  return(dat)
}

#' @title Fast Version of network::delete.vertices for Edgelist-formated Network
#'        Modified version for migrations.
#'
#' @description It will delete the edges, but will keep the same attribute number
#'    because when migration happens, only edges get removed.
#'
#' @param el A two-column matrix of current edges (edgelist) with an attribute
#'           variable \code{n} containing the total current network size.
#' @param vid A vector of IDs to delete from the edgelist.
#'
#' @details
#' This function is used in \code{EpiModel} modules to remove vertices (nodes)
#' from the edgelist object to account for exits from the population (e.g.,
#' deaths and out-migration)
#'
#' @return
#' Returns a updated edgelist object, \code{el}, with the edges of deleted
#' vertices removed from the edgelist and the ID numbers of the remaining edges
#' permuted downward.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("EpiModel")
#' set.seed(12345)
#' nw <- network_initialize(100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
#'
#' # Set seed for reproducibility
#' set.seed(123456)
#'
#' # networkLite representation structure after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#'
#' # Current edges
#' head(dat$el[[1]], 20)
#'
#' # Remove nodes 1 and 2
#' nodes.to.delete <- 1:2
#' dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes.to.delete)
#'
#' # Newly permuted edges
#' head(dat$el[[1]], 20)
#' }
#'
delete_edges <- function(el, vid) {

  new.el <- el
  if (length(vid) > 0) {
    el.rows.to.del <- which(el[, 1] %in% vid | el[, 2] %in% vid)
    if (length(el.rows.to.del) > 0) {
      new.el <- el[-el.rows.to.del, , drop = FALSE]
    }
    attributes(new.el)$n <- attributes(el)$n
  }

  return(new.el)
}


#' Save id and time of individuals that departure the network
#'
#' @description It aims to save the ID and time of infected individual that
#'  departure the network via natural or HIV related cause.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param departures ids of departures
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as departure_IDs.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
save_departures <- function(dat, prefix = NULL){

  if(!is.null(dat$stats$departures) == TRUE){
    if(is.null(prefix)){
      filename <- "departure_IDs.csv"
    } else {
      filename <- paste(prefix, "departure_IDs.csv", sep = "_")
    }

    write.csv(dat$stats$departures, file = filename, row.names = FALSE)
  }
}


#' Save id and time of individuals that departure the network
#'
#' @description It aims to save the ID and time of infected individual that
#'  departure the network via natural or HIV related cause.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param departures ids of departures
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as departure_IDs.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
set_departures <- function(dat, departures, at){
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  uid <- get_attr(dat, "unique_id")
  origin <- get_attr(dat, "origin")
  tx.status <- get_attr(dat, "tx.status")

  status_dep <- status[departures]
  stage_dep <- stage[departures]
  origin_dep <- origin[departures]
  tx.status_dep <- tx.status[departures]

  if(any(status_dep == "i")){
    #active_test <- active[departures]
    infID <- uid[departures]
    time <- at

    inf_time_df <- data.frame(time, infID, status_dep, stage_dep, origin_dep, tx.status_dep)
    #inf_time_df <- data.frame(time, infID, status_dep, stage_dep, origin_dep)
    inf_time_df <- subset(inf_time_df, status_dep == "i")

    #browser()
    if(!is.null(dat$stats$departures) == TRUE){
      inf_time_df <- rbind(dat$stats$departures, inf_time_df)
    }
    dat$stats$departures <- inf_time_df
  }

  return(dat)
}

#' Save stage of HIV infection of nodes in the network at final step
#'
#' @description At the final step of network simulation, it will save the IDs of
#'    infected nodes and their stage of HIV infection.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as stage_and_IDs.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
save_stage <- function(dat, prefix = NULL){


  #check whether the list stage_info exist
  if(!is.null(dat$stats$stage_info) == TRUE){
    if(is.null(prefix)){
      filename <- "stage_and_IDs.csv"
    } else {
      filename <- paste(prefix, "stage_and_IDs.csv", sep = "_")
    }

    write.csv(dat$stats$stage_info, file = filename, row.names = FALSE)
  }
}


#' Set stage of HIV infection of nodes in the network at final step
#'
#' @description At the final step of network simulation, it will save the IDs of
#'    infected nodes and their stage of HIV infection.
#'
#' @inheritParams EpiModel::arrivals.net
#'
#' @details
#' If a prefix is not provided, csv file will be saved as stage_and_IDs.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
set_stage <- function(dat, at){
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  diag.status <- get_attr(dat, "diag.status")
  uid <- get_attr(dat, "unique_id")
  tx.status <- get_attr(dat, "tx.status")

  #browser()
  status_inf_index <- which(status == "i" & active == 1)
  status_inf <- status[status_inf_index]
  stage_inf <- stage[status_inf_index]
  diag.status_inf <- diag.status[status_inf_index]
  tx.status_inf <- tx.status[status_inf_index]

  if(length(status_inf_index) > 0){
    active_test <- active[status_inf_index]
    infID <- uid[status_inf_index]

    inf_stage_df <- data.frame(infID, active_test, status_inf, stage_inf,
                               diag.status_inf, tx.status_inf)
    #inf_stage_df <- data.frame(infID, active_test, status_inf, stage_inf,
    #                           diag.status_inf)

    dat$stats$stage_info <- inf_stage_df
  }

  return(dat)

}



#' Save time in which antiretroviral treatment (ART) started
#'
#' This function will save the time in which an individual initiated ART
#'
#' @description Whenever an individual initiates ART, this function will save
#'    the IDs of infected nodes and when ART started.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as ART_init.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
save_art <- function(dat, prefix = NULL){

  if(!is.null(dat$stats$art_init) == TRUE){

    if(is.null(prefix)){
      filename <- "ART_init.csv"
    } else {
      filename <- paste(prefix, "ART_init.csv", sep = "_")
    }

    write.csv(dat$stats$art_init, file = filename, row.names = FALSE)
  }
}


#' Set time in which antiretroviral treatment (ART) started
#'
#' This function will save the time in which an individual initiated ART
#'
#' @description Whenever an individual initiates ART, this function will save
#'    the IDs of infected nodes and when ART started.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param IDs IDs of individuals who started ART.
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as ART.csv
#'
#' @return
#' @export
#'
set_art_init <- function(dat, at, IDs){
  #browser()
  art_init <- data.frame(time = at, IDs = IDs)

  if(!is.null(dat$stats$art_init) == TRUE){
    art_init <- rbind(dat$stats$art_init, art_init)
  }
  dat$stats$art_init <- art_init

  return(dat)
}


#' Save time in which antiretroviral treatment (ART) was hlated
#'
#' This function will save the time in which an individual initiated ART
#'
#' @description Whenever an individual initiates ART, this function will save
#'    the IDs of infected nodes and when ART started.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as ART_halt.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
save_art_halt <- function(dat, prefix = NULL){

  if(!is.null(dat$stats$art_halt) == TRUE){

    if(is.null(prefix)){
      filename <- "ART_halt.csv"
    } else {
      filename <- paste(prefix, "ART_halt.csv", sep = "_")
    }

    write.csv(dat$stats$art_halt, file = filename, row.names = FALSE)
  }
}


#' Set time in which antiretroviral treatment (ART) was halted
#'
#' This function will save the time in which an individual initiated ART
#'
#' @description Whenever an individual initiates ART, this function will save
#'    the IDs of infected nodes and when ART started.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param IDs IDs of individuals who started ART.
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as ART.csv
#'
#' @return
#' @export
#'
set_art_halt <- function(dat, at, IDs){
  #browser()
  art_halt <- data.frame(time = at, IDs = IDs)

  if(!is.null(dat$stats$art_halt) == TRUE){
    art_halt <- rbind(dat$stats$art_halt, art_halt)

  }
  dat$stats$art_halt <- art_halt

  return(dat)
}


#' Save time in which antiretroviral treatment (ART) was reinitiated
#'
#' This function will save the time in which an individual initiated ART
#'
#' @description Whenever an individual initiates ART, this function will save
#'    the IDs of infected nodes and when ART started.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as ART_reinit.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
save_art_reinit <- function(dat, prefix = NULL){

  if(!is.null(dat$stats$art_reinit) == TRUE){

    if(is.null(prefix)){
      filename <- "ART_reinit.csv"
    } else {
      filename <- paste(prefix, "ART_reinit.csv", sep = "_")
    }

    write.csv(dat$stats$art_reinit, file = filename, row.names = FALSE)
  }
}


#' Set time in which antiretroviral treatment (ART) was reinitiated
#'
#' This function will save the time in which an individual initiated ART
#'
#' @description Whenever an individual initiates ART, this function will save
#'    the IDs of infected nodes and when ART started.
#'
#' @inheritParams EpiModel::arrivals.net
#' @param IDs IDs of individuals who started ART.
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as ART.csv
#'
#' @return
#' @export
#'
set_art_reinit <- function(dat, at, IDs){
  #browser()
  art_reinit <- data.frame(time = at, IDs = IDs)

  if(!is.null(dat$stats$art_reinit) == TRUE){
    art_reinit <- rbind(dat$stats$art_reinit, art_reinit)
  }
  dat$stats$art_reinit <- art_reinit

  return(dat)
}


#'@export
set_transmat2 <- function (dat, del, at)
{
  del <- del[!duplicated(del$sus), ]
  del[["sus"]] <- get_unique_ids(dat, del[["sus"]])
  del[["inf"]] <- get_unique_ids(dat, del[["inf"]])
  browser()
  if (at != 2) {
    del <- rbind(dat$stats$transmat, del)
  }
  dat$stats$transmat <- del
  return(dat)
}

