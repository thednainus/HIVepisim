# Get transmission matrix ----
library(EpiModel)
library(HIVepisim)



#sim <- readRDS("Analyses/sim_for_tests.RDS")
sim <- readRDS("Analyses/Preliminary_results/results_tergmLite4/run_2/results_sim.RDS")
sim <- readRDS("Analyses/Preliminary_results/results_tergmLite1/run_2/results_sim.RDS")
sim <- readRDS("Analyses/Preliminary_results/results_tergmLite8/run_2/results_sim.RDS")

tm <- get_transmat(sim)
# to make sure that transmission only happens between individuals from same origin
tm[tm$infOrigin != tm$susOrigin,]
# get transmat by seed
tm2 <- get.transmat.phylo(tm, by_areas = "region", max_value = 501)


# testing codes for Matthew's program (VirusTreeSimulator) ----
# create transmission matrix file to be used with VirusTreeSimulator program
# if using the complete transmission matrix
create_inf_csv(tm, time_tr = 1, prefix="test")
# file with
create_sample_csv(tm, time_seq = 3652, seq_count = 1, prefix = "test")


# if using the list of transmission matrix by seed

 for(x in 1:length(tm2)){
   print(names(tm2[x]))
   print(paste(output_dir, names(tm2[x]), sep =""))

   create_inf_csv(tm2[[x]], time_tr = 1, prefix=paste(output_dir, names(tm2[x]), sep =""))
   create_sample_csv(tm2[[x]], time_seq = 3652, seq_count = 1, prefix = paste(output_dir, names(tm2[x]), sep =""))
 }

 # TO DO: Run bash script to run VirusTreeSimulator
