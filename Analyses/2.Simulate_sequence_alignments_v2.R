# Simulate sequence alignments
# Seq-Gen MUST be installed in your computer
library(stringr)
library(ape)

# Location for Seq-Gen. It should be changed to the correct location on your computer.
Software <- "/Applications/Seq-Gen-1.3.4/source/seq-gen"
#parameter for Seq-Gen
#seq-gen -mHKY -t0.5 -fe -l10000 -n1 < filename > filename+'.txt'
# values for generating sequence alignment for HIV pol as described in van der Kyl and Berkhout, 2012
# values of -s0.0028 based on paper by Patino-Galindo and Gonzalez-Candelas 2017
# values of -t8.75 based on Chen et al 2004
# Updated parameters values based on sequences for San Diego + background sequences


# get all VirusTreeSimulator in which branch lengths have been converted into years
# and tip names are in the format of ID_migrant
#vts_trees_years <- list.files(path = "output/vts", pattern ="*migrant_years", full.names = TRUE)
vts_trees_years <- list.files(path = "output/vts", pattern ="*onlyactive.tre", full.names = TRUE)

for(tree_name in vts_trees_years){
  filename_prefix <- str_split(tree_name, pattern = "\\.")[[1]][1]
  filename_prefix <- str_split(filename_prefix, pattern = "/")[[1]][3]

  #Create directory for sequence alignment
  if (!dir.exists("output/vts/alignments")) {
    dir.create("output/vts/alignments")
  }



  # simulate sequence alignment of 1000bp-
  seq_filename1000 <- paste("output/vts/alignments/", filename_prefix, "_ali_1000bp.fasta", sep = "")
  seq_filename10000 <- paste("output/vts/alignments/", filename_prefix, "_ali_10000bp.fasta", sep = "")


  # simulate sequence alignment using Seq-Gen
  # Run VirusTreeSimulator
  #cmd <- paste(parameters1, parameters2, vts_trees_years[1])
  #Simulate alignments of 1000 bp
  args1 <- c("-mHKY", "-l1000", "-t8.75", "-f0.389,0.165,0.228,0.218", "-s0.0028", "-n1", "-of", tree_name)
  system2(command = Software, args = args1, stdout = seq_filename1000)
  #Simulate alignments of 10000 bp
  args2 <- c("-mHKY", "-l10000", "-t8.75", "-f0.389,0.165,0.228,0.218", "-s0.0028", "-n1", "-of", tree_name)
  system2(command = Software, args = args2, stdout = seq_filename10000)
}
