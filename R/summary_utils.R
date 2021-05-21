#Summary utils: useful for plotting


#' Summarize data for statistics calculated for infector probability (W)
#'
#'
#' @param df dataframe containing data to summarize
#'
#' @return dataframe for median and quantiles for quantities of interest.
#'   see \code{link{validadeW}}
#' @export
#'
summarize_data <- function(df){

  df <- df[c(1,5,9,18,19)]

  #convert df to long format
  df_long <- melt(df, id.vars = c("sim", "n_tips_region"))

  #calculate median by simulation
  median_by_sim <- aggregate(x = df_long$value ,                # Specify data column
                             by = list(df_long$sim, df_long$n_tips_region, df_long$variable),       # Specify group indicator
                             FUN = median)
  names(median_by_sim)[4] <- "median"

  #lower quantile
  lowerq_by_sim <- aggregate(x = df_long$value ,
                                by = list(df_long$sim, df_long$n_tips_region, df_long$variable),
                                FUN = quantile,
                                probs=0.025)
  names(lowerq_by_sim)[4] <- "lower"


  #upper quantile
  upperq_by_sim <- aggregate(x = df_long$value ,
                             by = list(df_long$sim, df_long$n_tips_region, df_long$variable),
                             FUN = quantile,
                             probs=0.975)
  names(upperq_by_sim)[4] <- "upper"

  all_data <- cbind(median_by_sim, lowerq_by_sim[4], upperq_by_sim[4])

  return(all_data)

}



#' Get epi data for simulation runs
#'
#' @param list_dirs list of directories containing data
#'
#' @return A list of dataframes for each type of simulations.
#'
#' @details This function will go over each simulation run saved as RDS format,
#'   convert it to a dataframe and assign as simulation number the number of the
#'   respective run. This function was designed when running more than 1
#'   simulation independently in the cluster, to merge all the runs together for
#'   plotting purposes.
#'
#' @export
#'
get_epi_data <- function(list_dirs){

  #browser()

  for(n in 1:length(list_dirs)){

    if(n == 1){
      sim_all <- list()
    }

    type_name <- SplitPath(list_dirs[n])$filename
    list_runs <- dir(list_dirs[n], full.names = TRUE)

    count = 1
    for(r in 1:length(list_runs)){

      rds_sim <- list.files(list_runs[r], pattern = ".RDS", full.names = TRUE)
      #covert network sim into dataframe
      sim <- readRDS(rds_sim)
      #convert sim to data.frame
      sim_df <- as.data.frame(sim)
      sim_df$sim <- count
      sim_df$type <- type_name
      count <- count + 1

      if(r == 1){
        sim_dfs <- sim_df
      } else{sim_dfs <- rbind(sim_dfs, sim_df)}
    }
    sim_all[[n]] <- sim_dfs
  }
  return(sim_all)
}


#' Get rates to construct ROC curve
#'
#' Get false positive rates and true positive rate for ROC curves
#'
#' @param threshold list object with threshold values rangng from 0 to 1
#' @param df_true data.frame object for the infector probability and the true
#'    classification of transmission pairs
#'
#' @details The true data for df_true object will have 1 for a transmission pair
#'    that have occurred and 0 for a transmission pair that did not occur. This
#'    will be based by comparision with the transmission matrix.
#'
#' @return data.frame object for values of threshold, true positive rate (TPR)
#'    and false positive rate (FPR).
#' @export
get_rates <- function(threshold, df_true){

  df_true["pred"] <- ifelse(df_true$infectorProbability >= threshold, "1", "0")
  df_true$pred <- as.factor(df_true$pred)

  #create confusion matrix
  cm <- confusionMatrix(df_true$pred, df_true$labels, positive = "1")
  TPR <- cm$byClass[[1]]
  FPR <- 1 - cm$byClass[[2]]

  rates <- data.frame(threshold = threshold, TPR = TPR, FPR = FPR)

  return(rates)
}
