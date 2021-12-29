
#' @title Dynamic Network Updates
#'
#' @description This function handles all calls to the network object contained
#'              on the master dat object handled in \code{netsim}..
#'
#' @param dat Master list object containing a full \code{networkDynamic} object
#'        or networkLite edgelist (if using tergmLite), and other initialization
#'        information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#'
nwupdate_mig <- function(dat, at) {

  ## Attributes
  type <- get_control(dat, "type", override.null.error = TRUE)
  tergmLite <- get_control(dat, "tergmLite")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  active <- get_attr(dat, "active")
  entrTime <- get_attr(dat, "entrTime")
  exitTime <- get_attr(dat, "exitTime")
  migrationTime <- get_attr(dat, "migrationTime")
  migrants <- get_attr(dat, "migrant")
  origin <- get_attr(dat, "origin")

  final_step <- get_control(dat, "when2save_at")

  #statOnNw <- "status" %in% dat$temp$nwterms
  resimulate.network <- get_control(dat, "resimulate.network")
  isTERGM <- get_control(dat, "isTERGM")
  save.nwstats <- get_control(dat, "save.nwstats")
  save_nodes <- get_control(dat, "save_nodes")

  ## Vital Dynamics
  arrivals <- which(active == 1 & entrTime == at)
  departures <- which(active == 0 & exitTime == at)
  migrations12 <- which(migrants == 12 & active == 1 & migrationTime == at)
  migrations21 <- which(migrants == 21 & active == 1 & migrationTime == at)

  nArrivals <- length(arrivals)
  if (nArrivals > 0) {

    ## Arrivals--------
    nwterms <- dat$temp$nwterms
    if (!is.null(nwterms)) {
      curr.tab <- get_attr_prop(dat, nwterms)
      dat <- auto_update_attr(dat, arrivals, curr.tab)
    }
    if (length(unique(sapply(dat$attr, length))) != 1) {
      stop("Attribute list of unequal length. Check arrivals.net module.\n",
           print(cbind(sapply(get_attr_list(dat), length))))
    }
    if (tergmLite == FALSE) {
      dat$nw[[1]] <- add.vertices(dat$nw[[1]], nv = nArrivals)

      if (isTERGM == TRUE) {

        dat$nw[[1]] <- activate.vertices(dat$nw[[1]], onset = at,
                                         terminus = Inf, v = arrivals)
        dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                                 prefix = "testatus",
                                                 value = status[arrivals],
                                                 onset = at, terminus = Inf,
                                                 v = arrivals)
        dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]], prefix = "global_track",
                                                 value = origin[arrivals], onset = at,
                                                 terminus = Inf, v = arrivals)
      }
    }
    if (tergmLite == TRUE) {
      dat$el[[1]] <- add_vertices(dat$el[[1]], nv = nArrivals)
    }

  }


  ## Departures----------
  if (length(departures) > 0) {
    #if("i" %in% dat$attr$status[departures]){
      #browser()
    #}

    # save a csv file for time of departure and ID of infected individual that
    # is not in the network anymore
    # to be used with phylogenetics for getting time of terminal branch correctly.
    #save_departures(dat, departures, at)
    dat <- set_departures(dat, departures, at)

    if (tergmLite == FALSE) {
      if (isTERGM == TRUE) {
        dat$nw[[1]] <- deactivate.vertices(dat$nw[[1]], onset = at,
                                           terminus = Inf, v = departures,
                                           deactivate.edges = TRUE)
      } else {
        dat$nw[[1]] <- delete.vertices(dat$nw[[1]], vid = departures)
        dat <- delete_attr(dat, departures)
      }
    }

    if (tergmLite == TRUE) {
      dat <- delete_attr(dat, departures)
      dat$el[[1]] <- delete_vertices(dat$el[[1]], departures)

      if (dat$control$tergmLite.track.duration) {
        dat$p[[1]]$state$nw0 %n% "lasttoggle" <-
          delete_vertices(dat$p[[1]]$state$nw0 %n% "lasttoggle", departures)
      }
    }
  }

  ## Migrations---------
  if (length(migrations12) > 0) {
    #dat <- track_origin(dat, at, "12", migrations21)
    if (tergmLite == FALSE) {
        if(isTERGM == TRUE) {
          e = NULL
          for (vert in migrations12){
            e = c(e, get.edgeIDs.active(dat$nw[[1]], v = vert, onset = at,
                                        terminus = Inf, neighborhood = "combined"))
          }
          if (length(e) > 0) {
            dat$nw[[1]] <- deactivate.edges(dat$nw[[1]], onset = at, terminus = Inf,
                             e = unique(e))
          }
        }
    }

    if (tergmLite == TRUE) {
      dat$el[[1]] <- delete_edges(dat$el[[1]], migrations12)
    }
  }

  if (length(migrations21) > 0) {
    #dat <- track_origin(dat, at, "21", migrations21)

    if (tergmLite == FALSE) {
      if(isTERGM == TRUE) {
        e = NULL
        for (vert in migrations21){
          e = c(e, get.edgeIDs.active(x = dat$nw[[1]], v = vert, onset = at,
                                      terminus = Inf, neighborhood = "combined"))
        }
        if (length(e) > 0) {
          dat$nw[[1]] <- deactivate.edges(dat$nw[[1]], onset = at, terminus = Inf,
                                          e = unique(e))
        }

      }
    }

    if (tergmLite == TRUE) {
      dat$el[[1]] <- delete_edges(dat$el[[1]], migrations21)
    }
  }

  ## Update temporally extended disease and origin status
  if (tergmLite == FALSE & isTERGM == TRUE) {
    #idsNewInf <- which(status == "i" & infTime == at)
    #if (length(idsNewInf) > 0) {
    dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                             prefix = "testatus",
                                             value = status,
                                             onset = at,
                                             terminus = Inf)
    #}
  }

  # migrations
  if (tergmLite == FALSE & isTERGM == TRUE) {
    #idsNewMigs <- which((migrants == 12 | migrants == 21) & migrationTime == at)
    #if (length(idsNewMigs) > 0) {
    dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                             prefix = "global_track",
                                             value = origin,
                                             onset = at,
                                             terminus = Inf)
    #}
  }



  # save final list of infected individuals and origin
  #if(at == final_step){
  #  save_stage(dat)
  #}



  #browser()

  if (save.nwstats == TRUE){
    #get number of nodes by attribute origin
    total_nodes <- matrix(table(get.vertex.attribute(dat$nw[[1]], attrname = "origin")),
                          nrow = 1, ncol = 2)
    colnames(total_nodes) <- c("global", "region")
    dat <- track_nodes(dat, at, total_nodes)
  }

  # Save info for stage of HIV infection at the end of simulation
  # Save info for stage of HIV infection at the end of simulation
  if (at == final_step & dat$control$save.stats == TRUE){
    #browser()
    dat <- set_stage(dat, at)
    save_stage(dat, prefix = NULL)
    save_art(dat, prefix = NULL)
    save_art_halt(dat, prefix = NULL)
    save_art_reinit(dat, prefix = NULL)
    save_departures(dat, prefix = NULL)
    save_diagnosis_time(dat, prefix = NULL)
    save_track_stages(dat, prefix = NULL)
    save_track_origin(dat, prefix = NULL)
    save_r0(dat, prefix = NULL)

    if(save_nodes == TRUE){
      save_nodes(dat, prefix = NULL)
    }

  }



  ## Copy static attributes to network object
  if (tergmLite == FALSE & resimulate.network == TRUE) {
    dat <- copy_datattr_to_nwattr(dat)
  }


  # Record network in nw_list for x-sect ERGM simulations
  if (tergmLite == FALSE & isTERGM == FALSE & resimulate.network == TRUE) {
    dat$temp$nw_list[[at]] <- dat$nw[[1]]
  }
  if (tergmLite == FALSE & isTERGM == FALSE & resimulate.network == FALSE) {
    dat$temp$nw_list[[at]] <- set_vertex_attribute(dat$temp$nw_list[[at]],
                                                   "status", status)
  }

  # Attribute consistency checks
   if (tergmLite == FALSE) {
     tst <- get.vertex.attribute.active(dat$nw[[1]], "testatus", at = at)
     tst2 <- get.vertex.attribute.active(dat$nw[[1]], "global_track", at = at)
     if (any(is.na(tst))) {
       stop("Error in nwupdate.net: NA's in testatus attribute.\n")
     }
     if (any(is.na(tst2))) {
       stop("Error in nwupdate.net: NA's in global_track attribute.\n")
     }
     if (statOnNw == TRUE) {
       fstat <- get_vertex_attribute(dat$nw[[1]], "status")
       fstat2 <- get_vertex_attribute(dat$nw[[1]], "origin")
       if (!identical(tst, fstat)) {
         stop("Error in nwupdate.net: mismatch between status and testatus
               attribute.\n")
       }
       if (!identical(tst2, fstat2)) {
         stop("Error in nwupdate.net: mismatch between global_track and origin
               attribute.\n")
       }
     }
   }

  ## Output
  return(dat)
}
