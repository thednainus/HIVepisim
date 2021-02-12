#' @title Convert transmat infection tree into a phylo object
#'
#' @description Converts the edgelist matrix in the \code{transmat} object into
#'              a \code{phylo} object by doing the required reordering and
#'              labeling.
#'
#' @param x An object of class \code{"transmat"}, the output from
#'        \code{\link{get_transmat}}.
#' @param format If origin return tip in the form of ID_global or ID_region.
#'    If migrant return tip in the form of ID_1, ID_2, ID_21, ID_12
#' @param collapse.singles logical, DEPRECATED
#' @param vertex.exit.times  optional numeric vector providing the time of
#'        departure of vertices, to be used to scale the lengths of branches
#'        reaching to the tips. Index position on vector corresponds to network
#'        id. NA indicates no departure, so branch will extend to the end of the
#'        tree.
#' @param ...  further arguments (unused)
#'
#' @details
#' Converts a \code{\link{transmat}} object containing information about the
#' history of a simulated infection into a \code{\link{phylo}} object
#' representation suitable for plotting as a tree with
#' \code{\link[ape]{plot.phylo}}. Each infection event becomes a 'node'
#' (horizontal branch) in the resulting phylo tree, and each network vertex
#' becomes a 'tip' of the tree. The infection events are labeled with the vertex
#' id of the infector to make it possible to trace the path of infection.
#'
#' The infection timing information is included to position the phylo-nodes,
#' with the lines to the tips drawn to the max time value +1 (unless
#' \code{vertex.exit.times} are passed in it effectively assumes all vertices
#' are active/alive until the end of the simulation).
#'
#' If the transmat contains multiple infection seeds (there are multiple trees
#' with seperate root nodes) it will return a list of class 'multiPhylo', each
#' element of which is a phylo object.  See \code{\link[ape]{read.tree}}.
#'
#' Note that in EpiModel versions <= 1.2.4, the phylo tree was constructed
#' differently, translating network vertices to both phylo-nodes and tips and
#' requiring 'collapse.singles' to prune it to an appropriate branching
#' structure.
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.5)
#' init <- init.net(i.num = 1)
#' control <- control.net(type = "SI", nsteps = 40, nsims = 1, verbose = FALSE)
#'
#' mod1 <- netsim(est1, param, init, control)
#' tm <- get_transmat(mod1)
#' tmPhylo <- as.phylo.transmat(tm)
#' plot(tmPhylo, show.node.label = TRUE,
#'               root.edge = TRUE,
#'               cex = 0.5)
#'
toPhylo_transmatOrigin <- function(x, format = "migrant",
                              collapse.singles,
                              vertex.exit.times,
                              ...) {
  format <-  format
  # warnings if someone tries to use old args that are no longer supported
  if (!missing(collapse.singles)) {
    warning("the 'collapse.singles' argument to as.phylo.transmat is no longer
            supported and will be ignored")
  }

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
    sub_phylos <- lapply(v, function(v_sub) {
      # walk down the list to find elements below v_sub
      sub_rows <- which(tm$inf == v_sub)
      toFind <- v_sub
      while (length(toFind) > 0) {
        i <- toFind[1]
        sub_rows <- unique(c(sub_rows, which(tm$inf == i)))
        toFind <- c(toFind[-1], tm$sus[which(tm$inf == i)])
      }
      # call as.phylo on the subset of the edgelist
      toPhylo_transmatOrigin(tm[sub_rows, , drop = FALSE], format = format,
                        vertex.exit.times = vertex.exit.times)

    })
    names(sub_phylos) <- paste("seed", v, sep = "_")
    class(sub_phylos) <- c("multiPhylo", class(sub_phylos))
    return(sub_phylos)
  }

  el <- cbind(tm$inf, tm$sus)
  origNodes <- unique(as.vector(el))

  if (format == "origin"){
    el_ori <- cbind(at = tm$at,
                    inf = tm$inf,
                    sus = tm$sus,
                    inf_ori = paste(tm$inf, tm$infOrigin, sep = "_"),
                    sus_ori = paste(tm$sus, tm$susOrigin, sep = "_"))
  }

  if (format == "migrant"){
    el_ori <- cbind(at = tm$at,
                    inf = tm$inf,
                    sus = tm$sus,
                    inf_ori = paste(tm$inf, tm$infMigrant, sep = "_"),
                    sus_ori = paste(tm$sus, tm$susMigrant, sep = "_"))
  }


  # get tip names in the form of inf_infOri or sus_susOri
  tipNamesOri <- unique(as.vector(el_ori[,c(4:5)]))
  # get tip names only without the origin. This will allow to determine duplications
  # when a tip was region and changed to global (or vice versa, for example)
  tip_names <- unlist(lapply(tipNamesOri, function(x) str_split(string = x, pattern = "_")[[1]][1]))
  # get the duplicated
  dup_tip <- tip_names[duplicated(tip_names)]

  if(length(dup_tip) > 0){
    tip_origin_names <- get_newTip_names(dup_tip, el_ori)
    final_tip_names <- replace_dups(dup_tip, origNodes, tipNamesOri, tip_origin_names, tip_names)
  } else {
    final_tip_names <- tipNamesOri
  }


  # translate ids in el to sequential integers starting from one
  el[, 1] <- match(el[, 1], origNodes)
  el[, 2] <- match(el[, 2], origNodes)


  maxTip <- max(el)  # need to know what phylo node ids will start

  maxTime <- max(x$at) + 1
  if (!is.null(vertex.exit.times)) {
    maxTime <- max(maxTime, vertex.exit.times, na.rm = TRUE)
  }
  # create new ids for phyloNodes
  phyloNodes <- seq(from = maxTip + 1, length.out = length(origNodes) - 1)
  Nnode <- length(phyloNodes)
  # create labels for each phylo node based on infector id
  #phylo.label <- tm$inf
  # this is an alternate label form like i_j
  #phylo.label <- sapply(1:length(phyloNodes),function(r) {
  #  paste(tm[r,"inf"],tm[r,"infOrigin"],sep="_")
  #})
  phylo.label <- sapply(1:length(phyloNodes),function(r) {
    paste(tm[r,"inf"],tm[r,"sus"],sep="_")
  })

  # set default durations
  # since we don't know how long the graph vertices live, assume entire duration
  durations <- rep(NA, length(phyloNodes) * 2)
  tipExitTimes <- rep(maxTime, maxTip)
  if (!is.null(vertex.exit.times)) {
    # replace any NA values with max time
    vertex.exit.times[is.na(vertex.exit.times)] <- maxTime + 1
    # copy the vertex exit times into the appropriate positions in the
    # durations array
    durations[seq_len(maxTip)] <- vertex.exit.times[origNodes]
    # reorder the vertex.exit times to match new ids of tips
    tipExitTimes <- vertex.exit.times[origNodes]
  }

  # create a new edgelist by stepping through the existing edgelist
  # and creating the new links from phylo nodes to graph vertices (tips)
  # and from phylo node to phylo node
  # have to do this as progressive modifications

  # assume at least one xmit has occured
  # create the phylo node linking to the first
  # infector and infectee
  phyloEl <- rbind(cbind(phyloNodes[1], el[1, 1]),
                   cbind(phyloNodes[1], el[1, 2]))

  durations[1] <- tipExitTimes[el[1, 1]] - tm[["at"]][1]
  durations[2] <- tipExitTimes[el[1, 2]] - tm[["at"]][1]

  phyloN <- 1
  # loop over remaining rows
  if (nrow(el) > 1) {
    for (r in 2:nrow(el)) {
      # find id of infector
      infector <- el[r, 1]
      # find the phylo row of phylo node corresponding to the infector
      phyNRow <- which(phyloEl[, 2] == infector)
      # replace the infector with a new phylo node
      phyloEl[phyNRow, 2] <- phyloNodes[phyloN + 1]
      # link the new phylo node to the infector
      phyloEl <- rbind(phyloEl, cbind(phyloNodes[phyloN + 1], infector))
      # link the new phylo node to the infectee (tip)
      phyloEl <- rbind(phyloEl, cbind(phyloNodes[phyloN + 1], el[r, 2]))

      # update the timing on the replaced row that linked to tip
      durations[phyNRow] <-
        durations[phyNRow] - (tipExitTimes[infector] - tm[["at"]][r])
      # add timings for new rows equal to remaining time
      # infector
      durations[nrow(phyloEl) - 1] <- tipExitTimes[infector] - tm[["at"]][r]
      # infectee
      durations[nrow(phyloEl)] <- tipExitTimes[el[r, 2]] - tm[["at"]][r]


      # increment the phylo node counter
      phyloN <- phyloN + 1
    }
  }

  # format the output
  out <- list()
  out[["edge"]] <- phyloEl
  out[["Nnode"]] <- Nnode  # number of non-tip nodes
  out[["tip.label"]] <- final_tip_names
  out[["node.label"]] <- phylo.label
  out[["root.edge"]] <- x$at[1] # have to assume sim started at 0
  out[["edge.length"]] <- durations

  class(out) <- "phylo"
  return(out)
}



# function to get the oldest time
#' @export
oldest_time <- function(dup_rows, get_last_time){
  dup_rows <- as.data.frame(dup_rows)
  row_int <- dup_rows[dup_rows$at == get_last_time,]

  return(row_int)
}

# function to get the name of the tip in the for of id_origin
# because there can be duplicates, we need to remove those from the tree
#' @export
name_of_tips <- function(dup_tip, row_int){

  if(dup_tip %in% row_int$inf){
    infs_values <- row_int$inf_ori
  }else if(dup_tip %in% row_int$sus){
    infs_values <- row_int$sus_ori
  }

  return(infs_values)

}


#replace duplicated by the correct tip name
#' @export
replace_dups <- function(dup_tip, origNodes, tipNamesOri, tip_names_toInsert, tip_names){
  #get index of location of duplicated sequences
  index_to_remove <- which(tip_names %in% dup_tip)

  tip_tmp <- tipNamesOri

  tipNamesOri <- tip_tmp[-index_to_remove]
  # add new tip names
  #index to position to add tips
  index_to_add <- match(dup_tip, origNodes)
  index_to_add <- index_to_add - 1
  if(length(tip_names_toInsert) != length(index_to_add)){
    stop("length of `tip_names_toInsert` should be the same as `index_to_add`")
  }
  while(length(index_to_add) > 0){
    tipNamesOri2 <- append(x = tipNamesOri, values = tip_names_toInsert[1], after = index_to_add[1])
    tipNamesOri <- tipNamesOri2
    index_to_add <- index_to_add[-1]
    tip_names_toInsert <- tip_names_toInsert[-1]
  }

  return(tipNamesOri2)
}


#' @export
get_newTip_names <- function(dup_tip, el_ori){
  if(length(dup_tip) == 1){
    # get duplicated rows in the el_ori matrix
    dup_rows <- el_ori[el_ori[,2] == dup_tip | el_ori[,3] == dup_tip,]
    get_last_time <- as.numeric(tail(sort(dup_rows[,1]), 1))
    #row_int <- dup_rows[dup_rows$at == get_last_time,]
    row_int <- oldest_time(dup_rows, get_last_time)

    tip_origin_names <- name_of_tips(dup_tip, row_int)

  } else if (length(dup_tip) > 1 ){
    # get duplicated rows in the el_ori matrix
    dup_rows <- lapply(dup_tip, function(x) el_ori[el_ori[,2] == x | el_ori[,3] == x,])
    # if duplicated tip appears in more the one row, get the oldest time
    get_last_time <- lapply(dup_rows, function(x) as.numeric(tail(sort(x[,1]), 1)))
    #get the row of interest
    row_int <- mapply(oldest_time, dup_rows, get_last_time, SIMPLIFY = FALSE)

    tip_origin_names <- unname(mapply(name_of_tips, dup_tip, row_int))
  }

  return(tip_origin_names)

}
