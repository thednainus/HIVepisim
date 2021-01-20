#' @title Get transmat infection tree for each seed
#'
#' @description Separate a \code{transmat} object into different seedling
#'       event.
#'
#' @param x An object of class \code{"transmat"}, the output from
#'        \code{\link{get_transmat}}.
#' @param vertex.exit.times  optional numeric vector providing the time of
#'        departure of vertices, to be used to scale the lengths of branches
#'        reaching to the tips. Index position on vector corresponds to network
#'        id. NA indicates no departure, so branch will extend to the end of the
#'        tree.
#'
#' @details
#' Converts a \code{\link{transmat}} object containing information about the
#' history of a simulated infection into a separated transmission matrix by
#' seedling event.
#'
#' This function was created using as base the function
#' \code{\link{as.phylo.transmat}}
#'
#'
#' @export
#'

get.transmat.phylo <- function(x, vertex.exit.times) {


  # if not named properly, assume inf, sus at
  if (!all(c("inf", "sus", "at") %in% names(x))) {
    warning("input does not have appropriate column names for transmat,
            assuming first 3 should be 'inf','sus','at'")
    names(x) <- c("inf", "sus", "at")
  }
  tm <- x
  if (missing(vertex.exit.times)) {
    vertex.exit.times <- NULL
  }
  # find roots (infectors that never appear as sus)
  v <- setdiff(unique(tm$inf), unique(tm$sus))
  if (length(v) > 1) {
    message("found multiple trees, returning a list of ", length(v),
            "phylo objects")
    # need to extract the portions of the edgelist and call seperately
    subtm <- lapply(v, function(v_sub) {
      # walk down the list to find elements below v_sub
      sub_rows <- which(tm$inf == v_sub)
      toFind <- v_sub
      while (length(toFind) > 0) {
        i <- toFind[1]
        sub_rows <- unique(c(sub_rows, which(tm$inf == i)))
        toFind <- c(toFind[-1], tm$sus[which(tm$inf == i)])
      }
      # call as.phylo on the subset of the edgelist
      get.transmat.phylo(tm[sub_rows, , drop = FALSE],
                         vertex.exit.times = vertex.exit.times)
      if (exists("subtm")){
        subtm <- rbind(subtm, tm[sub_rows, , drop = FALSE])
      }else{
        subtm <- tm[sub_rows, , drop = FALSE]
      }

    })
    names(subtm) <- paste("seed", v, sep = "_")
    return(subtm)
  }
}


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
  dat <- update_uids(dat, n.new)
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
update_uids <- function(dat, n.new) {
  last_uid <- if (is.null(dat[["_last_uid"]])) 0L else dat[["_last_uid"]]
  next_uids <- seq_len(n.new) + last_uid
  dat[["_last_uid"]] <- last_uid + as.integer(n.new)
  dat <- append_attr(dat, "uid", next_uids, n.new)

  return(dat)
}



#' Save origin of nodes in the network at final step
#'
#' @description At the final step of network simulation, it will save the IDs of
#'    infected nodes and their origin. Origin can take the value of region or
#'    global.
#'
#' @inheritParams EpiModel::arrivals.net
#' @inheritParams create_sample_csv
#'
#' @details
#' If a prefix is not provided, csv file will be saved as infected_origin.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
save_origin <- function(dat, prefix = NULL){
  active <- get_attr(dat, "active")
  origin <- get_attr(dat, "origin")
  status <- get_attr(dat, "status")
  uid <- get_attr(dat, "uid")

  infID <- uid[active == 1 & status == "i"]
  infIDindex <- match(infID, uid)
  infOrigin <- origin[infIDindex]

  inf_df <- data.frame(infID, infOrigin)

  if(is.null(prefix)){
    filename <- "infected_origin.csv"
  } else {
    filename <- paste(prefix, "infected_origin.csv", sep = "_")
  }

  write.csv(inf_df, file = filename, row.names = FALSE)
}


#' Create transmission matrix csv
#'
#' @description Create a transmission matrix file to be used with the
#'    VirusTreeSimulator. VirusTreeSimulator requires that seeds are included
#'    in the file.
#'
#' @inheritParams create_sample_csv
#' @param time_tr Time of transmission for seeds. The other time of transmission
#'    should be provided in the tm dataframe.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as inf.csv
#'
#' @return
#' @export
#'
#' @examples
#' TO DO
create_inf_csv <- function(tm, time_tr, prefix = NULL){

  seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

  if(length(time_tr) == 1){
    seed_idtr <- data.frame(seed_names, rep(NA, length(seed_names)), rep(time_tr, length(seed_names)))
  }else if(length(time_tr) == length(seed_names)){
    seed_idtr <- data.frame(seed_names, rep(NA, length(seed_names)), time_tr)
  }else{
    stop("`time_tr` should be of length 1 or length of number of seeds")
  }


  colnames(seed_idtr) <- c("IDREC", "IDTR", "TIME_TR")

  inf_sus <- data.frame(tm$sus, tm$inf, tm$at)
  colnames(inf_sus) <- c("IDREC", "IDTR", "TIME_TR")

  all_data <- rbind(seed_idtr, inf_sus)

  if(is.null(prefix)){
    filename <- "inf.csv"
  } else {
    filename <- paste(prefix, "inf.csv", sep = "_")
  }

  write_csv(all_data, filename)
}





#' Get sample csv file
#'
#' @param tm Transmission matrix as returned using the function \link[EpiModel]{get_transmat}
#' @param time_seq Vector for sample time. If one value is provided, it will be replicated
#'    to the size of samples in the transmission matrix (tm).
#' @param seq_count The number of sequences per sample in the transmission matrix.
#'    If one value is provided, it will be replicated to the size of samples in
#'    the transmission matrix (tm).
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as sample.csv
#'
#' @return
#' @export
#'
#' @examples
#' To Do
create_sample_csv <- function(tm, time_seq, seq_count, prefix = NULL){

  IDPOP <- union(tm$inf, tm$sus)

  if(length(time_seq) == 1){
    TIME_SEQ <- rep(time_seq, length(IDPOP))
  } else if (length(time_seq) == length(IDPOP)){
    TIME_SEQ <- time_seq
  }else if (length(time_seq) != length(IDPOP)){
    stop("`time_seq` should be of length 1 or length of samples in the transmission matrix")
  }


  if(length(seq_count) == 1){
    SEQ_COUNT <- rep(seq_count, length(IDPOP))
  } else if (length(seq_count) == length(IDPOP)){
    SEQ_COUNT <- seq_count
  }else if(length(seq_count) != length(IDPOP)){
    stop(" `seq_count` should be of length 1 or length of samples in the transmission matrix")
  }


  all_data <- data.frame(IDPOP, TIME_SEQ, SEQ_COUNT)

  if(is.null(prefix)){
    filename <- "sample.csv"
  } else {
    filename <- paste(prefix, "sample.csv", sep = "_")
  }


  write_csv(all_data, filename)

}

