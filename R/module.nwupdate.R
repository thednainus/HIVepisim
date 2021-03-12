
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

  final_step <- get_control(dat, "nsteps")

  statOnNw <- "status" %in% dat$temp$nwterms
  resimulate.network <- get_control(dat, "resimulate.network")

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
    if (tergmLite == TRUE) {
      dat$el[[1]] <- add_vertices(dat$el[[1]], nv = nArrivals)
    }

  }


  ## Departures----------
  if (length(departures) > 0) {
    # save a csv file for time of departure and ID of infected individual that
    # is not in the network anymore
    # to be used with phylogenetics for getting time of terminal branch correctly.
    save_departures(dat, departures, at)
    if (tergmLite == FALSE) {
      dat$nw[[1]] <- deactivate.vertices(dat$nw[[1]], onset = at,
                                         terminus = Inf, v = departures,
                                         deactivate.edges = TRUE)
    }
    if (tergmLite == TRUE) {
      dat <- delete_attr(dat, departures)
      dat$el[[1]] <- delete_vertices(dat$el[[1]], departures)
    }
  }

  ## Migrations---------
  if (length(migrations12) > 0) {
    if (tergmLite == FALSE) {
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
    if (tergmLite == TRUE) {
      dat$el[[1]] <- delete_edges(dat$el[[1]], migrations12)
    }
  }

  if (length(migrations21) > 0) {
    if (tergmLite == FALSE) {
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
    if (tergmLite == TRUE) {
      dat$el[[1]] <- delete_edges(dat$el[[1]], migrations21)
    }
  }

  ## Infection-----------
  if (tergmLite == FALSE) {
    idsNewInf <- which(status == "i" & infTime == at)
    if (length(idsNewInf) > 0) {
      dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]], prefix = "testatus",
                                               value = "i", onset = at,
                                               terminus = Inf, v = idsNewInf)
    }
  }

  # migrations
  if (tergmLite == FALSE) {
    idsNewMigs <- which((migrants == 12 | migrants == 21) & migrationTime == at)
    if (length(idsNewMigs) > 0) {
      dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]], prefix = "global_track",
                                               value = origin[idsNewMigs], onset = at,
                                               terminus = Inf, v = idsNewMigs)
    }
  }



  # save final list of infected individuals and origin
  if(at == final_step){
    save_stage(dat)
  }


  ## Copy static attributes to network object
  if (tergmLite == FALSE & resimulate.network == TRUE) {
    dat <- copy_datattr_to_nwattr(dat)
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
