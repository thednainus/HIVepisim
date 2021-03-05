# Get transmission matrix ----
# and run VirusTreeSimulator
# This code has been tested on a Mac OS and might not work on windows or linux
library(EpiModel)
library(HIVepisim)
library(DescTools)
library(stringr)

# This function will generate input file to be used with program
# VirusTreeSimulator
# and will also run VirusTreeSimulator for each combination of
# inf and sample file
# It will create a directory "output" if directory does not exist
# and will save results to the directory "output"

# Location for VirusTreeSimulator should be changed accordingly
Software <- "java -jar /Applications/VirusTreeSimulator/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#parameter for VirusTreeSimulator
parameters <- "-demoModel Constant -N0 1"




get_tm <- function(list_dirs, years = 40, area = "region", max_value = 501){

  for(x in 1:length(list_dirs)){

    dir_names <- dir(list_dirs[x], pattern = "run*", full.names = TRUE)


    for (n in 1:length(dir_names)){
      pathname <- SplitPath(dir_names[n])$normpath
      output_name <- paste(pathname, "output", sep="/")

      #Create directory named output if it does not exist
      if (!dir.exists(output_name)) {
        dir.create(output_name)
      }

      rds_file <- list.files(dir_names[n], pattern = "*.RDS", full.names = TRUE)

      # get list of infected departures
      dep <- read.csv(list.files(dir_names[n], pattern = "*.csv", full.names = TRUE))

      sim <- readRDS(rds_file)
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
          #print(names(tm2[x]))
          output <- paste(output_name, names(tm2[x]), sep ="/")
          output <- paste(output, PLWHIV, newinf, sep = "_")



          #inf file name for VirusTreeSimulator
          inf_file <- paste(output, "_inf.csv", sep = "")
          #sample file name for VirusTreeSimulator
          sample_file <- paste(output, "_sample.csv", sep = "")

          create_inf_csv(tm2[[x]],
                        time_tr = 1,
                        prefix=output)

          IDPOP <- union(tm2[[x]]$inf, tm2[[x]]$sus)

          #reorder dep$infID to IDPO
          index <- match(IDPOP, dep$infID)
          time_seqs <- dep$time[index]
          time_seqs[is.na(time_seqs)] <- years * 365

          create_sample_csv(tm2[[x]],
                            time_seq = time_seqs,
                            seq_count = 1,
                            prefix = output)

          #Create directory named VTS (for VirusTreeSimulator) if it does not exist
          output_name_vts <- paste(pathname, "output/vts/", sep = "/")
          if (!dir.exists(output_name_vts)) {
            dir.create(output_name_vts)
          }

          #prefix with location of output directory to save results of VirusTreeSimulator
          prefix_vts <-  paste(output_name_vts, "results_vts", sep = "")

          # Run VirusTreeSimulator
          cmd <- paste(Software, parameters, inf_file, sample_file, prefix_vts, sep = " ")
          system(cmd)

          # change tip names to the format ID_migrant
          tip_names <- toPhylo_transmatOrigin(tm2[[x]], format = "migrant",
                                              by_areas = area, max_value = max_value,
                                              tips_only = TRUE)
          pattern = "_1_simple.nex"
          seed_name_tm <- strsplit(seed_names[x], split = "_")[[1]][2]
          seed_name_tm <- paste("ID", seed_name_tm, sep = "_")
          vts_tree <- read.nexus(paste(prefix_vts, seed_name_tm, pattern, sep=""))

          # run VirusTreeSimulator
        }
      }
    }
  }
}

# get dir names for directory ending with ".pbs"; example of directory name
# of interest "3124512[10].pbs"
# These are the directory generated when submitting array jobs in the
# Imperial College cluster.
list_dirs <- dir("Analyses/Preliminary_results", full.names = TRUE)
get_tm(list_dirs)


# Run VirusTreeSimulator on output imput files generated with this script








