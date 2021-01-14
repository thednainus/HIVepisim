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

get.transmat.phylo_mod <- function(x, vertex.exit.times) {


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

