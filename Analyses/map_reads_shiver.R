# map simulated Illumina reads using shiver (https://github.com/ChrisHIV/shiver)

library(DescTools)
#library(reticulate)

#set R to use python2
#use_python("/Users/user/opt/miniconda2/bin/python2")

#prefix for path names
prefix <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/HIVepisim/"

# list directories
iva_outputs <- dir(path = "output_deepseq/vts/merged_trees/Illumina_reads",
                  full.names = TRUE)
iva_outputs <- iva_outputs[1]

#location of shiver
shiver <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/shiver/"

# location for shiver align contigs
shiver_align_contigs <- paste(shiver, "shiver_align_contigs.sh", sep = "")

# location for shiver map reads
shiver_map_reads <- paste(shiver, "shiver_map_reads.sh", sep = "")

#location for MyInitDir
MyInitDir <- paste(shiver, "MyInitDir", sep = "")

#location of config.sh file
config_file <- paste(shiver, "config_MyChanges.sh", sep = "")





for(i in 1:length(iva_outputs)){
  #list ID directories
  id_dirs <- list.files(path = iva_outputs[i], full.names = TRUE)

  for(n in 1:length(id_dirs)){

    #make new dir to run it shiver analysis in its own directory
    shiver_output <- paste(id_dirs[n], "shiver", sep = "/")
    if (!dir.exists(shiver_output)) {
      dir.create(shiver_output)
    }

    shiver_output <- paste(prefix, shiver_output, sep = "")

    #run shiver_align_contigs
    SID_shiver <- SplitPath(id_dirs[n])$filename
    output_IVA_name <- paste(SID_shiver, "output", sep = "_")
    contigs_fasta_iva <- paste(paste(SplitPath(id_dirs[n])$normpath, output_IVA_name, sep = "/"),
                               "contigs.fasta", sep = "/")

    parameters1 <- paste(MyInitDir, config_file, contigs_fasta_iva, SID_shiver, sep = " ")
    command_shiver <- paste("cd", shiver_output, "&&", "bash", shiver_align_contigs, sep = " ")
    shiver_and_args <- paste(command_shiver, parameters1, sep = " ")

    system(shiver_and_args)

    #run shiver_map_reads.sh
    #check if the better alignment file exists
    SID_wRefs <- paste(SID_shiver, "cut_wRefs.fasta", sep = "_")
    SID_wRefs_fullPath <- paste(shiver_output, SID_wRefs, sep = "/")

    if(file.exists(SID_wRefs_fullPath) == FALSE){
      SID_wRefs <- paste(SID_shiver, "raw_wRefs.fasta", sep = "_")
      SID_wRefs_fullPath <- paste(shiver_output, SID_wRefs, sep = "/")
    }

    SID_blast <- paste(SID_shiver, "blast", sep = ".")
    SID_blast_fullPath <- paste(shiver_output, SID_blast, sep = "/")

    reads_1 <- paste(SID_shiver, "1.fq", sep = "_")
    reads_1_fullPath <- paste(SplitPath(id_dirs[n])$normpath, reads_1, sep = "/")

    reads_2 <- paste(SID_shiver, "2.fq", sep = "_")
    reads_2_fullPath <- paste(SplitPath(id_dirs[n])$normpath, reads_2, sep = "/")

    parameters2 <- paste(MyInitDir, config_file, contigs_fasta_iva, SID_shiver,
                        SID_blast_fullPath, SID_wRefs_fullPath,
                        reads_1_fullPath, reads_2_fullPath,  sep = " ")

    command_shiver_map <- paste("cd", shiver_output, "&&", "bash", shiver_map_reads, sep = " ")
    shiver_and_args2 <- paste(command_shiver_map, parameters2, sep = " ")

    system(shiver_and_args2)

  }
}

