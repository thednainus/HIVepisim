#' Title Validate infector probability calculation
#'
#' This function saves several values to compare W calculated on the true trees
#'    and the true transmission matrix. For the return values see below.
#'
#' @param tm Transmission matrix per seed
#' @param W_true Infector probability calculated on the true trees
#' @param seed_ID Seed ID name. Seed is the node ID that started the
#'    transmissions; it is the initial infected node in the simulations.
#' @param MH maximum height in which infector probabilities were calculated
#' @param tree object of class phylo
#' @param prefix prefix to save results. If prefix = NULL, results will be saved
#'   on file W_on_true.csv.
#'
#' @return a data frame of 1 row with the following values:
#'  1. Code = code indicating if W was calculated on true trees
#'  2. Maximum height = maximum height (MH) in which calculations were performed;
#'  3. Total tips = total number of tips in the phylogenetic tree;
#'  4. n_tips_region = total number of tips from region in the phylogenetic tree;
#'  5. n_tips_global = total number of tips from global in the phylogenetic tree;
#'  6. n_trans_W = total number of transmissions independent of infector probability
#'     values;
#'  7. n_trans_tm = total number of transmissions based on transmission matrix
#'     subset by MH;
#'  8. n_true_trans_W = total number of true transmission independent of infector
#'     probability (W) values. Here is the true transmission when compared to the
#'     transmission matrix;
#'  9. n_W80_all = number of transmission in which W >= 80%;
#'  10. n_W80_correctDonorRecipt = number of transmission in point 9 that is a
#'      true transmission (when comparing to the transmission matrix);
#'  11. n_W80_swapDonorRecipt = number of transmission in point 9 that is a
#'      true transmission (when comparing to the transmission matrix), but
#'      there was a swap between donor and recipient;
#'  12. expected_n_trans = expected number of transmission defined as the sum
#'      of all infector probability (W) values.
#' @export
#'
validadeW <- function(tm, W_true, W_estimated, seed_ID, MH, true_tree, estimated_tree, prefix = NULL){
  #get seed name
  IDnumber <- str_split(string = seed_ID, pattern = "_")[[1]][2]
  IDnumber <- paste("seed", IDnumber, sep = "_")

  #get tm by seed
  tm_seed <- tm[[IDnumber]]

  #add time in years
  tm_seed["time_years"] <- tm_seed["at"] * 1/365

  #subset tm to maximum height and to region only
  tm_mh <- subset(tm_seed, time_years >= MH &
                           infOrigin == "region" &
                           susOrigin == "region")

  # W on true trees ----
  W_on_true <- W_manipulations(W_true, code = "W on true")


  # W on estimated trees ----
  # converting W_estimated on a dataframe
  W_on_estimated <- W_manipulations(W_estimated, code = "W on estimated")


  #converting tm_mh to the same collumn names as Wsub
  tm <- data.frame(donor_ID = tm_mh$inf, recip_ID = tm_mh$sus,
                   infectorProbability = 1, Code = "True transmission")


  # get transmission
  all_trans_Wtrue <- semi_join(W_on_true$Wsub, tm, by = c("donor_ID", "recip_ID"))
  all_trans_West <- semi_join(W_on_estimated$Wsub, tm, by = c("donor_ID", "recip_ID"))




  if(!is.null(W_on_true$W80)){
    # true transmissions: correct identification of donor and recipient
    W80true_correctDonorRecip <- semi_join(W_on_true$W80, tm, by = c("donor_ID", "recip_ID"))
  }
  if(!is.null(W_on_true$W80_trunc)){
    # true transmission but incorrect identification of donor and recipient
    W80true_swapDonorRecip <- semi_join(W_on_true$W80_trunc, tm, by = c("donor_ID", "recip_ID"))
  }

  if(!is.null(W_on_estimated$W80)){
    # true transmissions: correct identification of donor and recipient
    # based on estimated phylogenetic trees
    W80est_correctDonorRecip <- semi_join(W_on_estimated$W80, tm, by = c("donor_ID", "recip_ID"))
  }
  if(!is.null(W_on_estimated$W80_trunc)){
    # true transmission but incorrect identification of donor and recipient
    # based on estimated phylogenetic trees
    W80est_swapDonorRecip <- semi_join(W_on_estimated$W80_trunc, tm, by = c("donor_ID", "recip_ID"))

  }



  #get names in phylogenetic trees
  region <- unlist(lapply(true_tree$tip.label, function(x) grepl("_1$", x) | grepl("_21", x)))
  global <- unlist(lapply(true_tree$tip.label, function(x) grepl("_2$", x) | grepl("_12", x)))


  all_data <- data.frame(maximum_height = MH,
                         total_tips = length(true_tree$tip.label),
                         n_tips_region = sum(region),
                         n_tips_global = sum(global),
                         n_trans_Wtrue = nrow(W_on_true$Wsub),
                         n_trans_West = nrow(W_on_estimated$Wsub),
                         n_trans_tm = nrow(tm),
                         n_true_trans_Wtrue = nrow(all_trans_Wtrue),
                         n_true_trans_West = nrow(all_trans_West),
                         n_W80true_all = ifelse(is.null(W_on_true$W80), 0, nrow(W_on_true$W80)),
                         n_W80true_correctDonorRecipt = ifelse(is.null(W_on_true$W80),
                                                               0, nrow(W80true_correctDonorRecip)),
                         n_W80true_swapDonorRecipt = ifelse(is.null(W_on_true$W80_trunc),
                                                            0, nrow(W80true_swapDonorRecip)),
                         n_W80est_all = ifelse(is.null(W_on_estimated$W80), 0, nrow(W_on_estimated$W80)),
                         n_W80est_correctDonorRecipt = ifelse(is.null(W_on_estimated$W80), 0,
                                                              nrow(W80est_correctDonorRecip)),
                         n_W80est_swapDonorRecipt = ifelse(is.null(W_on_estimated$W80_trunc), 0,
                                                           nrow(W80est_swapDonorRecip)),
                         expected_n_trans_Wtrue = sum(W_on_true$W$infectorProbability),
                         expected_n_trans_West = sum(W_on_estimated$W$infectorProbability))

  if(is.null(prefix)){
    filename <- "W_on_true.csv"
  } else {
    filename <- paste(prefix, "W_on_true.csv", sep = "_")
  }

  write.table(all_data, file = filename, append = TRUE, sep = ",",
              row.names = FALSE, col.names = !file.exists(filename))

}



W_manipulations <- function(W, code){

  #tm_mh$sus is the recipient
  #tm_mh$inf is the donor

  # W on true trees ----
  #convert W (infector probability matrix) to dataframe
  W <- as.data.frame(W)
  W["donor_ID"] <- unlist(lapply(W$donor, function(x) str_split(x, pattern = "_")[[1]][1]))
  W["recip_ID"] <- unlist(lapply(W$recip, function(x) str_split(x, pattern = "_")[[1]][1]))
  W["Code"] <- code
  W$donor_ID <- as.numeric(W$donor_ID)
  W$recip_ID <- as.numeric(W$recip_ID)

  # subseting data frame W to contain only donor_ID and recip_ID
  # this is comparable to the transmission matrix
  Wsub <- data.frame(donor_ID = W$donor_ID, recip_ID = W$recip_ID,
                          infectorProbability = W$infectorProbability, Code = W$Code)
  Wsub$donor_ID <- as.integer(Wsub$donor_ID)
  Wsub$recip_ID <- as.integer(Wsub$recip_ID)

  # get only transmissions in which infector probability is >= 80%
  W80 <- Wsub[Wsub$infectorProbability >= 0.8,]
  # modify the order of recipient and donor just to check whether identification of
  # donor was not correct
  W80_trunc <- data.frame(donor_ID = W80$recip_ID, recip_ID = W80$donor_ID)

  return(list(W = W, Wsub = Wsub, W80 = W80, W80_trunc = W80_trunc))


}





