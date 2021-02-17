# Get transmission matrix ----
library(EpiModel)
library(HIVepisim)



sim <- readRDS("Analyses/sim_for_tests.RDS")

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
 lapply(tm2, function(x) create_inf_csv(x, time_tr = 1, prefix=names(x)))

 for(x in 1:length(tm2)){
   print(names(tm2[x]))

   create_inf_csv(tm2[[x]], time_tr = 1, prefix=names(tm2[x]))
   create_sample_csv(tm2[[x]], time_seq = 3652, seq_count = 1, prefix = names(tm2[x]))
 }

 # TO DO: Run bash script to run VirusTreeSimulator
