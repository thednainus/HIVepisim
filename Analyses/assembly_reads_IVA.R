# assemble Illumina reads with IVA
# IVA: https://sanger-pathogens.github.io/iva/

library(DescTools)

library(reticulate)
# You must have a conda environment in you computer to run this script
# conda enviroment here was installed with Python 3.5
# and it has IVA installed as well
# numpy must also be installed in your conda enviroment for you to use
# this script
use_condaenv(condaenv = "iva-env5",
             required = TRUE)

# with py_config() one can check if everything is configured correctly
# and proceed with the R script
py_config()

# location of IVA within the conda environment
Software <- paste("/Users/user/opt/miniconda2/envs/iva-env5/bin/iva")

#number of threads to use with IVA
threads <- "-t 4"


# list directories
fastq_dirs <- dir(path = "output_deepseq/vts/merged_trees/Illumina_reads",
                  full.names = TRUE)
fastq_dirs <- fastq_dirs[1]

for(i in 1:length(fastq_dirs)){
  #list fasta files
  fastq_id_dirs <- list.files(path = fastq_dirs[i], full.names = TRUE)

  #save IVA results in same directory as fastq files

  for(n in 1:length(fastq_id_dirs)){

    outputname <- paste(SplitPath(fastq_id_dirs[n])$normpath,
                        paste(SplitPath(fastq_id_dirs[n])$filename, "output", sep = "_"),
                        sep = "/")

    reads_fwd <- paste(SplitPath(fastq_id_dirs[n])$normpath,
                       paste(SplitPath(fastq_id_dirs[n])$filename, "_1.fq", sep = ""),
                       sep = "/")
    reads_rev <- paste(SplitPath(fastq_id_dirs[n])$normpath,
                       paste(SplitPath(fastq_id_dirs[n])$filename, "_2.fq", sep = ""),
                       sep = "/")

    parameters <- c("-f", reads_fwd, "-r", reads_rev, outputname, threads)


    system2(command = "ulimit", args =  "-n 2048")
    system2(command = Software, args = parameters)




  }




}
