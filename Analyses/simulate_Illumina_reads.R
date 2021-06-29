# Simulate Illumina reads based on ART
# it requres ART installed in your computer
# https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm

library(DescTools)

# parameters to use ART
# -f the fold of read coverage to be simulated or number of reads/read pairs
#    generated for each amplicon
# -i the filename of input DNA/RNA reference
# -l the length of reads to be simulated
# -na do not output ALN alignment file
# -o the prefix of output filename
# -p indicate a paired-end read simulation or to generate reads from both ends
#    of amplicons
#    NOTE: art will automatically switch to a mate-pair simulation if the given
#    mean fragment size >= 2000
# -sam indicate to generate SAM alignment file
# -ss  The name of Illumina sequencing system of the built-in profile used for
#      simulation
# -ef  indicate to generate the zero sequencing errors SAM file as well the regular one
#    NOTE: the reads in the zero-error SAM file have the same alignment positions
#    as those in the regular SAM file, but have no sequencing errors


Software <- "/Applications/art_bin_MountRainier/art_illumina"




# list directories
fasta_dirs <- dir(path = "output_deepseq/vts/merged_trees/alignments/by_ID",
                  full.names = TRUE)

for(i in 1:length(fasta_dirs)){

  #list fasta files
  fasta_files <- list.files(path = fasta_dirs[i], pattern = ".fasta",
                            full.names = TRUE)

  fasta_files <- fasta_files[1:10]
  print(fasta_files)

  #get name of tree directory
  dir_name <- SplitPath(fasta_dirs[i])$filename

  # Create directory to save Illumina reads
  if (!dir.exists("output_deepseq/vts/merged_trees/Illumina_reads")) {
    dir.create("output_deepseq/vts/merged_trees/Illumina_reads")
  }

  where2save_reads <- paste("output_deepseq/vts/merged_trees/Illumina_reads",
                            dir_name, sep = "/")

  # Create directory for each type of phylogenetic tree
  # (if active + region for example)
  if (!dir.exists(where2save_reads)) {
    dir.create(where2save_reads)
  }

  for(f in 1:length(fasta_files)){

    # Create directory for save reads for each ID in the phylogenetic tree
    out_dir <- SplitPath(fasta_files[f])$filename
    if (!dir.exists(paste(where2save_reads, out_dir, sep = "/"))) {
      dir.create(paste(where2save_reads, out_dir, sep = "/"))
    }

    out_filename <- paste(SplitPath(fasta_files[f])$filename, "_", sep = "")
    input_filename <- paste("-i", fasta_files[f], sep = " ")
    output_filename <- paste("-o", paste(where2save_reads, out_dir, out_filename, sep = "/" ), sep= " ")
    parameters <- c("-ss MSv3", input_filename, output_filename, "-p", "-l 250",
                    "-m 500", "-s 48", "-f 10000", "-na")
    system2(command = Software, args = parameters)

  }


}
# TO DO
# now do assembly
# then use shiver
# then use phyloscanner





