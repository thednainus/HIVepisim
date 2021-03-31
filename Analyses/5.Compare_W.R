# Compare infector probabilities
library(phydynR)
library(treedater)
library(stringr)
library(ape)

#file name to load W
#type if W is calculated on TRUE or ESTIMATED tree
#tree to get number of tips
#simulation type: to keep track of transmission rate and migration rate

sumW <- function(file_name, type, tree, sim_type){
  load(file_name)
  taxa_size <- Ntip(tree)

  if(type == "true"){
    W <- W_true
    # sum of the infector probabilities for a given phylogenetic tree
    sumW <- sum(W$infectorProbability)
    # count the number of pairs in which an infector probability was calculated
    countW <- length(W$infectorProbability)

    results <- c(sim_type, type, taxa_size, sumW, countW)

  }

  if(type == "estimated"){
    W <- W_estimated
    # sum of the infector probabilities for a given phylogenetic tree
    sumW <- sum(W$infectorProbability)
    # count the number of pairs in which an infector probability was calculated
    countW <- length(W$infectorProbability)

    results <- c(sim_type, type, taxa_size, sumW, countW)

  }



}


#get true number of infector probabilities

load("Analyses/Preliminary_results/results_tergmLite1/run_1/output/vts/W/ID_1752_migrant_years_1_simple.RData")
load("Analyses/Preliminary_results/results_tergmLite1/run_1/output/vts/W/ID_413_migrant_years_1_simple.RData")

# list directories for each run
list_dirs <- dir("Analyses/Preliminary_results", full.names = TRUE)
#list_dirs <- dir("Analyses/Preliminary_results/10000pop/", full.names = TRUE)

W_onTrue <- paste("/output/vts/W")
W_onEstimated <- paste("/output/vts/W_estimated")

for(x in list_dirs){
  #list dirs for each run
  dir_names <- dir(list_dirs[x], pattern = "run*", full.names = TRUE)



  for(y in dir_names){

    W_onTrue_path <- paste(dir_names[y], W_onTrue, sep = "")
    W_files_onTrue <- list.files(W_onTrue_path, pattern = ".RData", full.names = TRUE)

    W_onEstimated_path <- paste(dir_names[y], W_onEstimated, sep = "")
    W_files_onEstimated <- list.files(W_onEstimated_path, pattern = ".RData", full.names = TRUE)

    if(length(W_files_onTrue) != length(W_files_onEstimated)){
      stop("Number of files for W calculated on the true tree should be the same
           as W calculated on the estimated ML and dated trees")
    }

    #load W estimated using true trees
    for(i in 1:length(W_files_onTrue)){
      # get ID number just to confirm I am analyzing files estimated with the same
      # simulated network
      get_filename <- strsplit(W_files_onTrue[1], split = "/")[[1]][8]
      get_ID_name <- str_extract(get_filename, pattern = "ID_\\d+_")

      #check if ID name on true tree is the same as in the estimated tree
      test_str <- str_detect(string = W_files_onEstimated[1], pattern = get_ID_name)
      if(test_str == FALSE){
        stop("Files for W ontrue tree are not in the same order as files for W
             on estimated tree")
      }

      #load W on true tree
      load(W_files_onTrue[1])
      load(W_files_onEstimated[1])





      print(W_files_onTrue[i])
      print(W_files_onEstimated[i])





    }



    }

  }



}



names(tm2)

tm_seed_1752 <- tm2$seed_1752
#convert time from days to years

tm_seed_1752["time_years"] <- tm_seed_1752$at * 1/365
#subset tm to maximum height
tm_seed_1752_MH <- subset(tm_seed_1752, time_years >= MH)


