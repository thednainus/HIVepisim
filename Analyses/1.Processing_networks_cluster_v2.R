# Get transmission matrix ----
# and run VirusTreeSimulator
# This code has been tested on a Mac OS and might not work on windows or linux
library(EpiModel)
library(HIVepisim)
library(DescTools)
library(stringr)
library(ape)
library(phydynR)
#library(adephylo)
library(castor)
library(dplyr)


#getwd()
# [1] "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/HIVepisim"
setwd("../../Preliminary_analyses/results_tergmLite1/run100_1/")

# This function will generate input file to be used with program
# VirusTreeSimulator
# and will also run VirusTreeSimulator for each combination of
# inf and sample file
# It will create a directory "output" if directory does not exist
# and will save results to the directory "output"

# Location for VirusTreeSimulator. It should be changed to the correct location on your computer.
#Software <- "java -jar ../../Programs/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
Software <- "java -jar /Applications/VirusTreeSimulator/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#parameter for VirusTreeSimulator
parameters <- "-demoModel Constant -N0 1"

years <-  40
area <-  "all"
max_value <-  NULL
# year of last sample dates in simulations
last_sample_date <- 2021.0
#maximum height
MH <- years

#Create directory named output if it does not exist
if (!dir.exists("output")) {
  dir.create("output")
}

# read departure ID files
dep <- read.csv("departure_IDs.csv")
#read stages for the other IDs
stages <- read.csv("stage_and_IDs.csv")

sim <- readRDS("results_sim.RDS")
sim_df <- as.data.frame(sim)
#get total number of people living with HIV at end of simulation
totalPLWHIV <- sum(sim_df$i.num.pop1[tail(sim_df$time, n=365)])/365
PLWHIV <- paste("PLWHIV", totalPLWHIV, sep="_")
#get number of new infections in the past year
newinf_per_year  <- sum(sim_df$incid.pop1[tail(sim_df$time, n=365)])
newinf <- paste("newinf", newinf_per_year, sep="_")

tm <- get_transmat(sim)

# get transmat by seed
# here I only get the ones that started in region
# later check all trees
#tm2 <- get.transmat.phylo(tm, by_areas = area, max_value = max_value)
#seed_names <- names(tm2)
# if using the list of transmission matrix by seed
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


      output <- paste("output", "output", sep ="/")
      output <- paste(output, PLWHIV, newinf, sep = "_")

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


      create_sample_csv(tm, time_seq = time_seqs, seq_count = 1, prefix = output)

      #Create directory named VTS (for VirusTreeSimulator) if it does not exist
      if (!dir.exists("output/vts/")) {
        dir.create("output/vts/")
      }

      #prefix with location of output directory to save results of VirusTreeSimulator
      prefix_vts <-  "output/vts/results_vts"

      # Run VirusTreeSimulator
      cmd <- paste(Software, parameters, inf_file, sample_file, prefix_vts, sep = " ")
      system(cmd)
    }

  #read all trees from vts
  list_trees <- dir("output/vts", pattern = "*_simple.nex", full.names = TRUE)
  trees <- lapply(list_trees, read.nexus)
  #add root.edge
  trees_rootedge <- lapply(trees, add_root_edge, total_sim_steps = 365 * years)
  vts_tree <- merge_trees(trees_rootedge)


  # convert branch lengths from days to years
  tree_years <- convert_branches(tree = vts_tree, scale = 1/365)

  #get tip names from VirusTreeSimulato tree
  tip_names_vts <- tree_years$tip.label

  # Change tip names from vts format to ID_migrant format
  # change the tips names returned by VirusTreeSimulator to those tip_names
  # in the form of ID_migrant
  tree_years$tip.label <- reorder_tip_names(tip_names, tip_names_vts)


  # save tree to simulate sequence alignment using Python script
  tree_filename <- paste(prefix_vts, "_merged_trees_allnodes.tre", sep="")
  write.tree(phy = tree_years, file = tree_filename)

  #drop tips in which node has died before end of simulation
  all_IDs <- as.numeric(unlist(lapply(tree_years$tip.label, function(x) str_split(x, "_")[[1]][1])))
  all_IDs <- all_IDs %in% dep$infID
  all_IDs <- setNames(all_IDs, tree_years$tip.label)
  toDrop <- all_IDs[all_IDs == TRUE]
  only_active_tree <- drop.tip(tree_years, names(toDrop))

  # save tree to simulate sequence alignment using Python script
  tree_filename <- paste(prefix_vts, "_merged_trees_onlyactive.tre", sep="")
  write.tree(phy = only_active_tree, file = tree_filename)


  # Remove sequences from "global" -----
  # to calculate infector probability
  migrant_ID <- unlist(lapply(only_active_tree$tip.label,
                              function(x) ifelse(str_split(x, "_")[[1]][2] == "1" | str_split(x, "_")[[1]][2] == "21",
                                                 FALSE, TRUE)))
  migrant_ID <- setNames(migrant_ID, only_active_tree$tip.label)
  toDrop <- migrant_ID[migrant_ID == TRUE]
  region_only_tree <- drop.tip(only_active_tree, names(toDrop))

  # drop the ones that have not been diagnosed
  migrant_IDs <- as.numeric(unlist(lapply(tree_years$tip.label, function(x) str_split(x, "_")[[1]][1])))
  active_diagStatus <- subset(stages, diag.status_inf == 0)
  migrant_IDs <- migrant_IDs %in% active_diagStatus$infID
  migrant_IDs <- setNames(migrant_IDs, tree_years$tip.label)
  toDrop <- migrant_IDs[migrant_IDs == TRUE]
  only_active_tree_diag <- drop.tip(region_only_tree, names(toDrop))




  # Calculate infector probability ----
  all_cd4s <- get_cd4s(IDPOP = IDPOP, df_actives = stages, df_departures = dep)
  names(all_cd4s) <- tip_names

  #get ehi (early HIV infection)
  #named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  ehis <- ifelse(all_cd4s == 1e3, TRUE, FALSE)


  #match IDs from phylogenetic tree to the departure.csv file
  index <- match(IDPOP, dep$infID)
  time_seqs <- dep$time[index]
  time_seqs[is.na(time_seqs)] <- years * 365

  sampleTimes <- round(time_seqs/365, digits = 3)
  sampleTimes <- years - sampleTimes
  # if sampleTimes == 0, it means that ID was active at the last time step
  # of network simulation
  # I assume that the most recent sample data is in object last_sample_date
  sampleTimes <- ifelse(sampleTimes != 0, last_sample_date - sampleTimes, last_sample_date)
  names(sampleTimes) <- tip_names


  # to calculate infector probabilities
  # sampleTimes: must use years
  # cd4s: named numeric vector, cd4 at time of sampling
  # ehi: named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  # numberPeopleLivingWithHIV: scalar
  # numberNewInfectionsPerYear: scalar

  W <- phylo.source.attribution.hiv.msm( only_active_tree, sampleTimes[only_active_tree$tip.label],
                                         cd4s = all_cd4s[only_active_tree$tip.label],
                                         ehi = ehis[only_active_tree$tip.label],
                                         numberPeopleLivingWithHIV  = totalPLWHIV,
                                         numberNewInfectionsPerYear = newinf_per_year,
                                         maxHeight = years,
                                         res = 1e3,
                                         treeErrorTol = Inf)

  W1 <- phylo.source.attribution.hiv.msm( only_active_tree_diag, sampleTimes[only_active_tree_diag$tip.label],
                                         cd4s = all_cd4s[only_active_tree_diag$tip.label],
                                         ehi = ehis[only_active_tree_diag$tip.label],
                                         numberPeopleLivingWithHIV  = totalPLWHIV,
                                         numberNewInfectionsPerYear = newinf_per_year,
                                         maxHeight = years,
                                         res = 1e3,
                                         treeErrorTol = Inf)

  #newtree resampled tree
  n <- as.integer(length(only_active_tree$tip.label) * 0.5)
  newtree <- keep.tip( only_active_tree, sample( only_active_tree$tip.label,
                                                 size = n , replace=FALSE))

  # save tree to simulate sequence alignment using Python script
  newtree_filename <- paste(prefix_vts, "_region_resampled.tre", sep="")
  write.tree(phy = newtree, file = newtree_filename)


  W2 <- phylo.source.attribution.hiv.msm( newtree, sampleTimes[newtree$tip.label],
                                         cd4s = all_cd4s[newtree$tip.label],
                                         ehi = ehis[newtree$tip.label],
                                         numberPeopleLivingWithHIV  = totalPLWHIV,
                                         numberNewInfectionsPerYear = newinf_per_year,
                                         maxHeight = years,
                                         res = 1e3,
                                         treeErrorTol = Inf)


  #newtree resampled tree for diagnosed+region only tree
  n1 <- as.integer(length(only_active_tree_diag$tip.label) * 0.1)
  newtree1 <- keep.tip( only_active_tree_diag, sample( only_active_tree_diag$tip.label,
                                                 size = n1 , replace=FALSE))

  # save tree to simulate sequence alignment using Python script
  newtree_filename <- paste(prefix_vts, "_region_resampled_diag.tre", sep="")
  write.tree(phy = newtree1, file = newtree_filename)


  W3 <- phylo.source.attribution.hiv.msm( newtree1, sampleTimes[newtree1$tip.label],
                                          cd4s = all_cd4s[newtree1$tip.label],
                                          ehi = ehis[newtree1$tip.label],
                                          numberPeopleLivingWithHIV  = totalPLWHIV,
                                          numberNewInfectionsPerYear = newinf_per_year,
                                          maxHeight = years,
                                          res = 1e3,
                                          treeErrorTol = Inf)



  #Create directory named W (to save everything related to infector probability)
  # if it does not exist
  if (!dir.exists("output/vts/W")) {
    dir.create("output/vts/W")
  }

  W_filename <- paste("output/vts/W/", "merged_trees", "_migrant_years_1_simple", ".RData", sep="")
  save(years, max_value, last_sample_date, tm, tree_years, region_only_tree, region_only_tree_diag,
       sampleTimes,all_cd4s, ehis, newinf_per_year, totalPLWHIV, W,W1, W2, W3, newtree, newtree1,
       file = W_filename)

  summaryW2(sim = "1", run = "1",
            tm = tm, W, ID = "all_tips", MH = years,
            tree = tree_years, prefix = paste(prefix, "all_region", sep = "_"), labels = TRUE)

  summaryW2(sim = "1", run = "1",
            tm = tm, W1, ID = "all_tips_diag", MH = years,
            tree = tree_years, prefix = paste(prefix, "all_region_diag", sep = "_"), labels = TRUE)

  summaryW2(sim = "1", run = "1",
            tm = tm, W2, ID = "subset_tips", MH = years,
            tree = newtree, prefix = paste(prefix, "region_subset_all", sep = "_"), labels = TRUE)

  summaryW2(sim = "1", run = "1",
            tm = tm, W3, ID = "subset_tips", MH = years,
            tree = newtree, prefix = paste(prefix, "region_subset_all_diag", sep = "_"), labels = TRUE)

}

