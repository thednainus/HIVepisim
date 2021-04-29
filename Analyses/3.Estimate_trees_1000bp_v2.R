# Estimate phylogenetic tree using IQ-TREE
# and convert branch lengths to unit of calendar time
library(treedater)
library(DescTools)
library(phydynR)
library(ape)


# You have to download IQ-TREE to run this script
# Change to the correct path of IQ-TREE on your computer
Software <- "/Applications/iqtree-2.1.2-MacOSX/bin/iqtree2"
maxCPU <- 3


#setwd("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/HIVepisim/Analyses/Preliminary_results/results_tergmLite1/run_9")


#list files
list_files <- list.files("output/vts/alignments", pattern = "_1000bp.fasta", full.names = TRUE)
# read alignments of 1000 bp

list_sampleTimes <- list.files("output/vts/W", pattern = ".RData", full.names = TRUE)

#make dir to save iqtree
if (!dir.exists("output/vts/alignments/iqtree_results_1000bp")) {
  dir.create("output/vts/alignments/iqtree_results_1000bp")
}

for(ali in list_files){


  ali_name <- paste("-s", ali, sep = " ")
  # -czb Collapse near zero branches, so that the final tree may be
  # multifurcating.
  iqtree_param <- c(ali_name, "-m HKY", "-T AUTO", "-ntmax", maxCPU, "-czb")
  #run iqtree
  system2(command = Software, args = iqtree_param)

  #get sample times
  #ali_ID_name <- str_extract(ali, pattern = "vtsID_\\d+_")
  #ali_ID_name <- str_extract(ali, pattern = "ID_\\d+_")

  #get sampleTimes equivalent to the tree estimated
  #get index
  #index <- which(str_detect(list_sampleTimes, pattern = ali_ID_name))

  # get sampleTimes, cd4s, ehis to estimate treedater and infector probabilities
  load(list_sampleTimes)

  #estimate treedated tree
  #read ML tree
  mltree <- read.tree(paste(ali, ".treefile", sep = ""))

  #check if tree has polytomies
  #if tree has polytomies resolve polytomis randomly
  if(is.binary(mltree) == FALSE){
    mltree <- multi2di(mltree)
    if(is.rooted(mltree) == TRUE){
      mltree <- unroot(mltree)
    }
  }




  #if(length(sampleTimes) != length(mltree$tip.label)){
  #  stop("length of sampleTimes and tip_label should be the same")

  #}

  #run treedater
  dated_tree <- dater(tre = mltree, sts = sampleTimes[mltree$tip.label], s = 1000)


  #drop tips of the tree that are from "global"
  # to calculate infector probability
  migrant_ID <- unlist(lapply(dated_tree$tip.label, function(x) ifelse(str_split(x, "_")[[1]][2] == 1 | str_split(x, "_")[[1]][2] == 21, FALSE, TRUE)))
  migrant_ID <- setNames(migrant_ID, dated_tree$tip.label)
  toDrop <- migrant_ID[migrant_ID == TRUE]
  region_only_dated_tree <- drop.tip(dated_tree, names(toDrop))

  # calculate infector probability (W) on the estimated tree
  # using same cd4s, ehis and MH as in the true trees
  # sampleTimes: must use years
  # cd4s: named numeric vector, cd4 at time of sampling
  # ehi: named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  # numberPeopleLivingWithHIV: scalar
  # numberNewInfectionsPerYear: scalar

  W_estimated <- phylo.source.attribution.hiv.msm( region_only_dated_tree, sampleTimes[region_only_dated_tree$tip.label],
                                                   cd4s = all_cd4s[region_only_dated_tree$tip.label],
                                                   ehi = ehis[region_only_dated_tree$tip.label],
                                                   numberPeopleLivingWithHIV  = totalPLWHIV,
                                                   numberNewInfectionsPerYear = newinf_per_year,
                                                   maxHeight = MH,
                                                   res = 1e3,
                                                   treeErrorTol = Inf)

  #Create directory named W_estimated (to save everything related to infector probability)
  # if it does not exist
  if (!dir.exists("output/vts/W_estimated")) {
    dir.create("output/vts/W_estimated")
  }

  W_filename <- paste("output/vts/W_estimated/", ali_ID_name, "migrant_years_1_simple_estimated", ".RData", sep="")
  save(years, max_value, last_sample_date, tm, region_only_dated_tree,
       sampleTimes, all_cd4s, ehis, newinf_per_year, totalPLWHIV, W_estimated,
       file = W_filename)
}

