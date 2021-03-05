#combine data frames
#to check plots for number of interest (like number of infected individuals)

combine_df <- function(list_dirs){

  for(x in 1:length(list_dirs)){

    dir_names <- dir(list_dirs[x], pattern = "run*", full.names = TRUE)

    pathname <- SplitPath(dir_names[x])$fullpath
    output_name <- paste(pathname, "output_df", sep="")

    #Create directory named output if it does not exist
    if (!dir.exists(output_name)) {
      dir.create(output_name)
    }

    df_list <- list()
    count_sim <- 1
    for (n in 1:length(dir_names)){
      rds_file <- list.files(dir_names[n], pattern = "*.RDS", full.names = TRUE)

      sim <- readRDS(rds_file)
      sim_df <- as.data.frame(sim)
      sim_df["sim"] <- count_sim

      df_list[[n]] <- sim_df
      count_sim <- count_sim + 1
    }
    df_final <- do.call(rbind, df_list)

    return(df_final)
  }
}

list_dirs <- dir("Analyses/Preliminary_results", full.names = TRUE)
sim_dfs <- combine_df(list_dirs)

median.sim <- aggregate(sim_dfs[, 3:ncol(sim_dfs)], list(sim_dfs$time), median)
median.sim.l <- melt(median.sim, id.vars="Group.1")
colnames(median.sim.l)[3] <- "Median"

lower.sim <- aggregate(sim_dfs[, 3:ncol(sim_dfs)], list(sim_dfs$time), quantile, probs=0.025, na.rm = TRUE)
lower.sim.l <- melt(lower.sim, id.vars="Group.1")
colnames(lower.sim.l)[3] <- "Lower"


upper.sim <- aggregate(sim_dfs[, 3:ncol(sim_dfs)], list(sim_dfs$time), quantile, probs=0.975, na.rm = TRUE)
upper.sim.l <- melt(upper.sim, id.vars="Group.1")
colnames(upper.sim.l)[3] <- "Upper"

all.data_df <- data.frame(median = median.sim.l$Median, lower = lower.sim.l$Lower, upper = upper.sim.l$Upper)
all.data_df <- cbind(all.data_df, time = median.sim.l$Group.1, variable = median.sim.l$variable)


ggplot(subset(all.data_df, variable == "s.num.pop2"), aes(x=time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.20) +
  geom_line(aes(y = median)) +
  xlab("Time (days)") +
  theme_bw()


ggplot(all.data_df, aes(x=time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.20) +
  geom_line(aes(y = median)) +
  facet_wrap( ~ variable, scales = "free")
  xlab("Time (days)") +
  theme_bw()
