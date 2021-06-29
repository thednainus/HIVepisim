# Get transmission matrix ----
# and run VirusTreeSimulator
# This code has been tested on a Mac OS and might not work on windows or linux
library(EpiModel)
library(HIVepisim)
library(DescTools)
library(stringr)
library(ape)
library(phydynR)
library(castor)
library(dplyr)


# This function will generate input file to be used with program
# VirusTreeSimulator
# and will also run VirusTreeSimulator for each combination of
# inf and sample file
# It will create a directory "output" if directory does not exist
# and will save results to the directory "output"

# Location for VirusTreeSimulator. It should be changed to the correct location on your computer.
#Software <- "java -jar VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
Software <- "java -jar /Applications/VirusTreeSimulator/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#parameter for VirusTreeSimulator
parameters <- "-demoModel Constant -N0 1"

years <-  40
area <-  "all"
max_value <-  NULL
# year of last sample dates in simulations
last_sample_date <- 2021.0


#Create directory named output_deepseq if it does not exist
if (!dir.exists("output_deepseq")) {
  dir.create("output_deepseq")
}

# read departure ID files
dep <- read.csv("departure_IDs.csv")
#read stages for the other IDs
stages <- read.csv("stage_and_IDs.csv")

sim <- readRDS("results_sim.RDS")
sim_df <- as.data.frame(sim)

tm <- get_transmat(sim)


if(!is.null(tm)){

    # Get tip names in the form of ID_migrant
    tip_names <- get_tip_names(tm, format = "migrant",
                                   by_areas = area, max_value = max_value,
                                   tips_only = TRUE)

    # check number of individual within "region"
    # region is code as 1 and 21
    tip_names_migrant_ID <- unlist(lapply(tip_names, function(x) str_split(x, "_")[[1]][2]))
    #total_tips_region <- sum(tip_names_migrant_ID == "1" | tip_names_migrant_ID == "21")

    if(length(tip_names_migrant_ID) > 0){


      output <- paste("output_deepseq", "output", sep ="/")
      #output <- paste(output, PLWHIV, newinf, sep = "_")

      #inf file name for VirusTreeSimulator
      inf_file <- paste(output, "_inf.csv", sep = "")
      #sample file name for VirusTreeSimulator
      sample_file <- paste(output, "_sample.csv", sep = "")

      seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

      create_inf_csv(tm, time_tr = rep(0, length(seed_names)), prefix=output)


      el <- cbind(tm$inf, tm$sus)
      IDPOP <- unique(as.vector(el))


      #match IDs from phylogenetic tree to the departure.csv file
      index <- match(IDPOP, dep$infID)
      time_seqs <- dep$time[index]
      time_seqs[is.na(time_seqs)] <- years * 365

      # if tx.status = 0 then sample 10 different viruses
      # if tx.status = 1 then sample 2 different viruses
      # note that this is just to start the simulations
      # but I should try to find real values for those

      dep["seq_count"] <- ifelse(dep$tx.status_dep == 0, 10, 2)
      stages["seq_count"] <- ifelse(stages$tx.status_inf == 0, 10, 2)

      index_active <- match(IDPOP, stages$infID)

      count_dep <- dep$seq_count[index]
      count_active <- stages$seq_count[index_active]
      count_total <- count_dep
      count_total[is.na(count_total)] <- count_active[!is.na(count_active)]

      #commented lines are just to check if code above is correct
      #test_dep_id <- dep$infID[index]
      #test_active_id <- stages$infID[index_active]
      #all_ids <- test_dep_id
      #all_ids[is.na(all_ids)] <- test_active_id[!is.na(test_active_id)]

      create_sample_csv(tm, time_seq = time_seqs, seq_count = count_total, prefix = output)

      #Create directory named VTS (for VirusTreeSimulator) if it does not exist
      if (!dir.exists("output_deepseq/vts/")) {
        dir.create("output_deepseq/vts/")
      }

      #prefix with location of output directory to save results of VirusTreeSimulator
      prefix_vts <-  "output_deepseq/vts/results_vts"

      # Run VirusTreeSimulator
      cmd <- paste(Software, parameters, inf_file, sample_file, prefix_vts, sep = " ")
      system(cmd)
    }

  #Create directory to save phylogenetic trees
  if (!dir.exists("output_deepseq/vts/merged_trees/")) {
    dir.create("output_deepseq/vts/merged_trees")
  }

  #prefix with location of output directory to save merged phylogenetic trees
  prefix_trees <-  "output_deepseq/vts/merged_trees/results"

  #read all trees from vts
  list_trees <- dir("output_deepseq/vts", pattern = "*_simple.nex", full.names = TRUE)
  trees <- lapply(list_trees, read.nexus)
  #add root.edge
  trees_rootedge <- lapply(trees, add_root_edge, total_sim_steps = 365 * years, root.edge_value = 0)
  vts_tree <- merge_trees(trees_rootedge)


  # convert branch lengths from days to years
  tree_years <- convert_branches(tree = vts_tree, scale = 1/365)

  #get tip names from VirusTreeSimulator tree
  tip_names_vts <- tree_years$tip.label

  # Change tip names from vts format to ID_migrant format
  # change the tips names returned by VirusTreeSimulator to those tip_names
  # in the form of ID_migrant
  # tip_names1 will have appended to the ID name if it is a migrant or not
  # (1 is from region; 2 is from global;
  # 12 migrated from region to global;
  # 21 migrated from global to region)
  tip_names1 <- reorder_tip_names(tip_names, tip_names_vts)

  #get sample number for each sequence
  sample_numbers <- unlist(lapply(tip_names_vts, function(x) paste("sample",
                                                                   str_split(x, pattern = "_")[[1]][4],
                                                                   sep = "_")
    ))


  #tree_years$tip.label <- paste(tip_names1, sample_numbers, "count_1", sep = "_")
  tree_years$tip.label <- paste("ID", tip_names1, sample_numbers, sep = "_")

  # save trees
  tree_filename <- paste(prefix_trees, "_merged_trees_allnodes_years.tre", sep="")
  write.tree(phy = tree_years, file = tree_filename)

  #drop tips in which node has died before end of simulation
  all_IDs <- as.numeric(unlist(lapply(tree_years$tip.label, function(x) str_split(x, "_")[[1]][2])))
  all_IDs <- all_IDs %in% dep$infID
  all_IDs <- setNames(all_IDs, tree_years$tip.label)
  toDrop <- all_IDs[all_IDs == TRUE]
  only_active_tree <- drop.tip(tree_years, names(toDrop))

  # save tree to simulate sequence alignment
  tree_filename <- paste(prefix_trees, "_merged_trees_onlyactive.tre", sep="")
  write.tree(phy = only_active_tree, file = tree_filename)


  # drop the IDs of individuals that have not been diagnosed
  IDs <- as.numeric(unlist(lapply(only_active_tree$tip.label, function(x) str_split(x, "_")[[1]][2])))
  active_diagStatus <- subset(stages, diag.status_inf == 0)
  IDs <- IDs %in% active_diagStatus$infID
  IDs <- setNames(IDs, only_active_tree$tip.label)
  toDrop <- IDs[IDs == TRUE]
  only_active_tree_diag <- drop.tip(only_active_tree, names(toDrop))

  # save tree to simulate sequence alignment
  tree_filename <- paste(prefix_trees, "_merged_trees_onlyactive_diag.tre", sep="")
  write.tree(phy = only_active_tree_diag, file = tree_filename)


  # Remove sequences from "global" -----
  # to calculate infector probability
  migrant_ID <- unlist(lapply(only_active_tree_diag$tip.label,
                              function(x) ifelse(str_split(x, "_")[[1]][3] == "1" | str_split(x, "_")[[1]][3] == "21",
                                                 FALSE, TRUE)))
  migrant_ID <- setNames(migrant_ID, only_active_tree_diag$tip.label)
  toDrop <- migrant_ID[migrant_ID == TRUE]
  region_only_tree <- drop.tip(only_active_tree_diag, names(toDrop))

  # save tree to simulate sequence alignment
  tree_filename <- paste(prefix_trees, "_merged_trees_active_region_diag.tre", sep="")
  write.tree(phy = region_only_tree, file = tree_filename)


}

