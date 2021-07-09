# script to run phyloscanner: https://github.com/BDI-pathogens/phyloscanner
# this script run the second part of phyloscanner:
# analyse trees generated from part 1

# phyloscanner dependencies
# https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh

library(DescTools)
library(reticulate)

#location of phyloscanner make analyse_trees R program
phyloscanner <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner"
analyse_trees <- paste(phyloscanner, "phyloscanner_analyse_trees.R", sep = "/")

# phyloscanner analyse_trees parameters that does not depend on location
# of results of running run_phyloscanner_part1.R

#label to name output results
output_label <- "results_"

# outgroup name
outgroup <- "--outgroupName C.BW.00.00BW07621.AF443088"

# multifurcation threshold as g = "guess"
multifurcation <- "--multifurcationThreshold g"

#normRefFileName: using
#update this later using the 2019 reference genome sequences I used with shiver
#check manual on how to create this normalization file
norm <- "--normRefFileName /Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner/InfoAndInputs/HIV_DistanceNormalisationOverGenome.csv"




#

# list directories
output_dirs <- dir(path = "output_deepseq/vts/merged_trees/Illumina_reads",
                   full.names = TRUE)
output_dirs <- output_dirs[1]

for(i in 1:length(output_dirs)){

  #phyloscanner directory
  pyloscanner_dir_fullPath <- paste(SplitPath(output_dirs[i])$normpath,
                                    "phyloscanner", sep = "/")



}
