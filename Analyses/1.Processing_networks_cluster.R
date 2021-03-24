# Get transmission matrix ----
# and run VirusTreeSimulator
# This code has been tested on a Mac OS and might not work on windows or linux
library(EpiModel)
library(HIVepisim)
library(DescTools)
library(stringr)
library(ape)
library(phydynR)

#getwd()
# [1] "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/HIVepisim"
#setwd("Analyses/Preliminary_results/results_tergmLite4/run_8/")

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
last_sample_date <- 2021
#maximum height
MH <- 20

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
totalPLWHIV <- sim_df$i.num.pop1[tail(sim_df$time, n=1)]
PLWHIV <- paste("PLWHIV", totalPLWHIV, sep="_")
#get number of new infections per year
newinf_per_year  <- sum(sim_df$incid.pop1)/years
newinf <- paste("newinf", newinf_per_year, sep="_")

tm <- get_transmat(sim)
# get transmat by seed
# here I only get the ones that started in region
# later check all trees
tm2 <- get.transmat.phylo(tm, by_areas = area, max_value = max_value)
seed_names <- names(tm2)
# if using the list of transmission matrix by seed
if(!is.null(tm2)){
  for(x in 1:length(tm2)){

    # Get tip names in the form of ID_migrant
    tip_names <- toPhylo_transmatOrigin(tm2[[x]], format = "migrant",
                                        by_areas = area, max_value = max_value,
                                        tips_only = TRUE)

    # check number of individual within "region"
    # region is code as 1
    tip_names_migrant_ID <- unlist(lapply(tip_names, function(x) str_split(x, "_")[[1]][2]))
    total_tips_region <- sum(tip_names_migrant_ID == "1")

    if(total_tips_region >= 100){

      #print(names(tm2[x]))
      output <- paste("output", names(tm2[x]), sep ="/")
      output <- paste(output, PLWHIV, newinf, sep = "_")

      #inf file name for VirusTreeSimulator
      inf_file <- paste(output, "_inf.csv", sep = "")
      #sample file name for VirusTreeSimulator
      sample_file <- paste(output, "_sample.csv", sep = "")

      create_inf_csv(tm2[[x]], time_tr = 1, prefix=output)


      el <- cbind(tm2[[x]]$inf, tm2[[x]]$sus)
      IDPOP <- unique(as.vector(el))


      #match IDs from phylogenetic tree to the departure.csv file
      index <- match(IDPOP, dep$infID)
      time_seqs <- dep$time[index]
      time_seqs[is.na(time_seqs)] <- years * 365

      create_sample_csv(tm2[[x]], time_seq = time_seqs, seq_count = 1, prefix = output)

      #Create directory named VTS (for VirusTreeSimulator) if it does not exist
      if (!dir.exists("output/vts/")) {
        dir.create("output/vts/")
      }

      #prefix with location of output directory to save results of VirusTreeSimulator
      prefix_vts <-  "output/vts/results_vts"

      # Run VirusTreeSimulator
      cmd <- paste(Software, parameters, inf_file, sample_file, prefix_vts, sep = " ")
      system(cmd)


      pattern = "_1_simple.nex"
      seed_name_tm <- strsplit(seed_names[x], split = "_")[[1]][2]
      seed_name_tm <- paste("ID", seed_name_tm, sep = "_")
      vts_tree <- read.nexus(paste(prefix_vts, seed_name_tm, pattern, sep=""))


      # convert branch lengths from days to years
      tree_years <- convert_branches(tree = vts_tree, scale = 1/365)

      #get tip names from VirusTreeSimulato tree
      tip_names_vts <- tree_years$tip.label

      # Change tip names from vts format to ID_migrant format
      # change the tips names returned by VirusTreeSimulator to those tip_names
      # in the form of ID_migrant
      tree_years$tip.label <- reorder_tip_names(tip_names, tip_names_vts)

      # save tree to simulate sequence alignment using Python script
      tree_filename <- paste(prefix_vts, seed_name_tm, "_migrant_years_1_simple", ".tre", sep="")
      write.tree(phy = tree_years, file = tree_filename)


      # Remove sequences from "global" -----
      # to calculate infector probability

      migrant_ID <- ifelse(unlist(lapply(tree_years$tip.label, function(x) str_split(x, "_")[[1]][2])) != 1, TRUE, FALSE)
      migrant_ID <- setNames(migrant_ID, tree_years$tip.label)
      toDrop <- migrant_ID[migrant_ID == TRUE]
      region_only_tree <- drop.tip(tree_years, names(toDrop))



      # Calculate infector probability ----
      all_cd4s <- get_cd4s(IDPOP = IDPOP, df_actives = stages, df_departures = dep)
      names(all_cd4s) <- tip_names

      #get ehi (early HIV infection)
      #named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
      ehis <- ifelse(all_cd4s == 1e3, TRUE, FALSE)

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

      W <- phylo.source.attribution.hiv.msm( region_only_tree, sampleTimes[region_only_tree$tip.label],
                                             cd4s = all_cd4s[region_only_tree$tip.label],
                                             ehi = ehis[region_only_tree$tip.label],
                                             numberPeopleLivingWithHIV  = totalPLWHIV,
                                             numberNewInfectionsPerYear = newinf_per_year,
                                             maxHeight = MH,
                                             res = 1e3,
                                             treeErrorTol = Inf)

      #Create directory named W (to save everything related to infector probability)
      # if it does not exist
      if (!dir.exists("output/vts/W")) {
        dir.create("output/vts/W")
      }

      W_filename <- paste("output/vts/W/", seed_name_tm, "_migrant_years_1_simple", ".RData", sep="")
      save(years, MH, max_value, last_sample_date, tm, tm2, tree_years,
           all_cd4s, ehis, newinf_per_year, totalPLWHIV, W,
           file = W_filename)

    }
  }
}

