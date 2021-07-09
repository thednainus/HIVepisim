# script to run phyloscanner: https://github.com/BDI-pathogens/phyloscanner
# this script run the first part of phyloscanner:
# construct trees from reads

# phyloscanner dependencies
# https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh

library(DescTools)
library(reticulate)

use_python("/Users/user/opt/miniconda2/bin/python2")
py_config()

#location of phyloscanner make tree python code
phyloscanner <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner"
make_trees <- paste(phyloscanner, "phyloscanner_make_trees.py", sep = "/")

# window parameters
window_width <- "150,"
window_overlap <- "52,"
start_pos <- "520,"
end_pos <- "9480"

#make tree parameters
auto_window_param <- paste("--auto-window-params", window_width, sep = " ")
auto_window_param <- paste(auto_window_param, window_overlap,
                           start_pos, end_pos, sep = "")
alignment_other_refs <- paste("--alignment-of-other-refs",
                              paste(phyloscanner, "InfoAndInputs",
                                    "2refs_HXB2_C.BW.fasta", sep = "/"), sep = " ")

paiwise_align_to <- paste("--pairwise-align-to",
                          "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
                          sep = " ")

merge_pairs <- "--merge-paired-reads"

#quality_trim_ends <- paste("--quality-trim-ends", "25", sep = " ")

#min_internal_quality <- paste("--min-internal-quality", "25", sep = " ")

excision_ref <- paste("--excision-ref", "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
                      sep = " " )

excision_coords <- paste("--excision-coords", "$(cat",
                         paste(phyloscanner, "InfoAndInputs",
                               "DrugResistancePositionsInHXB2.txt", sep = "/"),
                         ")", sep = " ")

merging_threshold <- paste("--merging-threshold-a", "1", sep = " ")

min_read_count <- paste("--min-read-count", "2", sep = " ")

#raxml option
x_raxml <- paste("--x-raxml", '"/Applications/standard-RAxML/raxmlHPC-PTHREADS-AVX -m GTRCAT -p 1 --no-seq-check"',
               sep = " ")






#Phyloscanner location

# list directories
output_dirs <- dir(path = "output_deepseq/vts/merged_trees/Illumina_reads",
                   full.names = TRUE)
output_dirs <- output_dirs[1]

for(i in 1:length(output_dirs)){

  # phyloscanner directory
  pyloscanner_dir_fullPath <- paste(SplitPath(output_dirs[i])$normpath,
                                    "phyloscanner", sep = "/")

  BamsRefsAndIds <- paste(pyloscanner_dir_fullPath, "BamsRefsAndIDs.csv", sep = "/")


  #parameters <- paste(BamsRefsAndIds, auto_window_param, alignment_other_refs,
  #              paiwise_align_to, merge_pairs, quality_trim_ends,
  #              min_internal_quality, excision_ref, excision_coords,
  #              merging_threshold, min_read_count, raxml, sep = " ")

  parameters <- paste(BamsRefsAndIds, auto_window_param, alignment_other_refs,
                      paiwise_align_to, merge_pairs,
                      excision_ref, excision_coords,
                      merging_threshold, min_read_count, x_raxml, sep = " ")


  command_make_trees <- paste("cd", pyloscanner_dir_fullPath, "&&", make_trees, sep = " ")
  makeTrees_and_args <- paste(command_make_trees, parameters, sep = " ")


  system(makeTrees_and_args)



}
