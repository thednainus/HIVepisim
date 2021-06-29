# Simulate sequence alignments
# Seq-Gen MUST be installed in your computer
library(stringr)
library(ape)

# Location for Seq-Gen. It should be changed to the correct location on your computer.
Software <- "/Applications/Seq-Gen-1.3.4/source/seq-gen"
#reference genomic sequence to give as input to Seg-Gen
ref_genome <- "~/Desktop/Imperial/newHIVproject-01Aug2020/reference_subtypeB/genome/reference_subtypeB_completeGenome.phy"
# file containing the number of trees in which to simulate alignments
number_trees <- "~/Desktop/Imperial/newHIVproject-01Aug2020/reference_subtypeB/genome/tree_number.txt"
# parameter for Seq-Gen
# seq-gen -mHKY -t0.5 -fe -n1 -k 1 < filename > filename+'.txt'
# -k parameter is to simulate sequences based on ancestral sequences.
# I will use HXB2 sequence to simulate HIV genomic sequences
# values for generating sequence alignment for HIV pol as described in van der Kyl and Berkhout, 2012
# values of -s0.0028 based on paper by Patino-Galindo and Gonzalez-Candelas 2017
# values of -t8.75 based on Chen et al 2004
# Updated parameters values based on sequences for San Diego + background sequences


# get all VirusTreeSimulator in which branch lengths have been converted into years
# and tip names are in the format of ID_migrant
vts_trees_years <- list.files(path = "output_deepseq/vts/merged_trees",
                              pattern ="*onlyactive.tre|*onlyactive_diag.tre|*active_region_diag.tre", full.names = TRUE)

for(tree_name in vts_trees_years){
  filename_prefix <- str_split(tree_name, pattern = "\\.")[[1]][1]
  filename_prefix <- str_split(filename_prefix, pattern = "/")[[1]][4]


  # merge files with reference genome sequence and tree sequence for output
  # to Seq_gen
  #Create directory for sequence alignment
  if (!dir.exists("output_deepseq/vts/merged_trees/input_SeqGen")) {
    dir.create("output_deepseq/vts/merged_trees/input_SeqGen")
  }

  input_seqgen <- paste("output_deepseq/vts/merged_trees/input_SeqGen",
                        "tmp_inputSeqGen.txt", sep = "/")

  system2(command = "cat", args = c(ref_genome, number_trees, tree_name),
          stdout = input_seqgen)


  #Create directory for sequence alignment
  if (!dir.exists("output_deepseq/vts/merged_trees/alignments")) {
    dir.create("output_deepseq/vts/merged_trees/alignments")
  }



  # simulate sequence alignment of 1000bp-
  seq_filename_genome <- paste("output_deepseq/vts/merged_trees/alignments/",
                               filename_prefix, "_ali_genome.fasta", sep = "")


  # simulate sequence alignment using Seq-Gen
  # Run VirusTreeSimulator
  #cmd <- paste(parameters1, parameters2, vts_trees_years[1])
  #Simulate complete HIV genomes
  args1 <- c("-mHKY", "-t8.75", "-f0.389,0.165,0.228,0.218", "-s0.0028", "-n1", "-of", "-k 1", input_seqgen)
  system2(command = Software, args = args1, stdout = seq_filename_genome)

  # read alignment and split each combination of sequences per individual
  # to simulate Illumina reads
  genome_seqs <- read.FASTA(seq_filename_genome)

  #get ID names
  ID_names <- unlist(lapply(names(genome_seqs), function(x) str_split(string = x, pattern = "_")[[1]][2]))

  #create dataframe with sequence name in phylogenetic tree and ID names
  seq_ID_names <- data.frame(seq_name = names(genome_seqs), ID_names = ID_names)
  seq_ID_names$ID_names <- as.factor(seq_ID_names$ID_names)
  seq_ID_split <- split(seq_ID_names, f = seq_ID_names$ID_names)

  #match id_names to sequence name in DNAbin
  index <- lapply(seq_ID_split, function(x) match(unname(unlist(x[[1]])), names(genome_seqs)))

  seq_by_ID <- lapply(index, function(x) genome_seqs[head(x, n=1):tail(x, n=1)])

  #save separate alignment in fasta format

  if (!dir.exists("output_deepseq/vts/merged_trees/alignments/by_ID")) {
    dir.create("output_deepseq/vts/merged_trees/alignments/by_ID")
  }

  dir_name_by_tree <- paste("output_deepseq/vts/merged_trees/alignments/by_ID", filename_prefix, sep = "/")

  if (!dir.exists(dir_name_by_tree)) {
    dir.create(dir_name_by_tree)
  }


  lapply(seq_by_ID, function(x) write.FASTA(x,
                                            file = paste(dir_name_by_tree, "/",
                                                         paste("ID_", str_split(names(x)[1], pattern = "_")[[1]][2], ".fasta", sep = ""), sep = "" )))

}
