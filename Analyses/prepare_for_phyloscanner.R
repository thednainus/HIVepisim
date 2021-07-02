# prepare files to run phyloscanner

# shiver output files were run in individual directories per ID (shiver should
# be run in individual directory per ID)

# in this script, I will move all files that are required to run phyloscanner to
# a separate directory.

library(DescTools)

#prefix for path names
prefix <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/HIVepisim/"

# list directories
output_dirs <- dir(path = "output_deepseq/vts/merged_trees/Illumina_reads",
                   full.names = TRUE)
output_dirs <- output_dirs[1]


# Move files from ID/shiver to a phyloscanner directory ----


for(i in 1:length(output_dirs)){

  # make phyloscanner directory if it still does not exist
  pyloscanner_dir <- paste(output_dirs[i], "phyloscanner", sep = "/")
  pyloscanner_dir_fullPath <- paste(prefix, pyloscanner_dir, sep = "")
  if (!dir.exists(pyloscanner_dir)) {
    dir.create(pyloscanner_dir)
  }


  #list ID directories
  id_dirs <- list.files(path = output_dirs[i], full.names = TRUE)

  for(n in 1:length(id_dirs)){

    #make new dir to run it shiver analysis in its own directory
    shiver_output <- paste(id_dirs[n], "shiver", sep = "/")
    shiver_output <- paste(prefix, shiver_output, sep = "")

    #run shiver_align_contigs
    SID_shiver <- SplitPath(id_dirs[n])$filename

    if(SID_shiver != "phyloscanner"){

      # copy files to run phyloscanner to new directory
      # specifically I will copy files that ends in ID_remap.bam or ID.bam (this
      # latter will be copied if _remap.bam does not exist). Similarly, I will
      # copy ID_remap_ref.fasta or ID_ref.fasta (the latter will be copied if
      # _remap option does not exist)

      #mv bam files
      bam_file <- paste(SID_shiver, "remap.bam", sep = "_")
      bam_file_fullPath <- paste(shiver_output, bam_file, sep = "/")

      # check if file exists
      if(file.exists(bam_file_fullPath) == FALSE){
        bam_file <- paste(SID_shiver, "bam", sep = ".")
        bam_file_fullPath <- paste(shiver_output, bam_file, sep = "/")
      }

      system(paste("mv", bam_file_fullPath, pyloscanner_dir_fullPath, sep = " "))


      #mv bam.bai files
      bam_bai_file <- paste(SID_shiver, "remap.bam.bai", sep = "_")
      bam_bai_file_fullPath <- paste(shiver_output, bam_bai_file, sep = "/")

      # check if file exists
      if(file.exists(bam_bai_file_fullPath) == FALSE){
        bam_bai_file <- paste(SID_shiver, "bam.bai", sep = ".")
        bam_bai_file_fullPath <- paste(shiver_output, bam_bai_file, sep = "/")
      }

      system(paste("mv", bam_bai_file_fullPath, pyloscanner_dir_fullPath, sep = " "))


      #mv ref.fasta files
      ref_fasta_file <- paste(SID_shiver, "remap_ref.fasta", sep = "_")
      ref_fasta_file_fullPath <- paste(shiver_output, ref_fasta_file, sep = "/")

      # check if file exists
      if(file.exists(ref_fasta_file_fullPath) == FALSE){
        ref_fasta_file <- paste(SID_shiver, "ref.fasta", sep = ".")
        ref_fasta_file_fullPath <- paste(shiver_output, ref_fasta_file, sep = "/")
      }

      system(paste("mv", ref_fasta_file_fullPath, pyloscanner_dir_fullPath, sep = " "))
    }
  }
}


# Create input csv file for phyloscanner ----

bam_files <- list.files(path = pyloscanner_dir_fullPath, pattern = "*.bam$")
ref_files <- list.files(path = pyloscanner_dir_fullPath, pattern = "*remap_ref.fasta")

if(length(bam_files) != length(ref_files)){
  stop("length of bam_files has to be equal to length of ref_files")
}

# create dataframe of bam and ref files

input_files <- data.frame(bam = bam_files, ref_fasta = ref_files)

# check that ID names in bam matches ID names in ref_fasta

input_files["bam_ids"] <- unlist(lapply(input_files$bam, function(x) str_split(x, pattern = "_")[[1]][2]))
input_files["ref_ids"] <- unlist(lapply(input_files$ref_fasta, function(x) str_split(x, pattern = "_")[[1]][2]))

# check that all IDs in order are equal to each other
if(any(input_files$bam_ids == input_files$ref_ids) == FALSE){
  stop("All or some IDs between bam and ref.fasta files do not match")
}

input_files["IDs"] <- paste("ID", input_files$bam_ids, sep = "_")

#paste full file name path to bam and ref.fasta files
input_files["bam"] <- paste(pyloscanner_dir_fullPath, input_files$bam, sep = "/")
input_files["ref_fasta"] <- paste(pyloscanner_dir_fullPath, input_files$ref_fasta, sep = "/")

#final dataframe
input_files <- input_files[c(1,2,5)]

#save csv file that will be used as input file for phyloscanner
filename <- paste(pyloscanner_dir_fullPath, "BamsRefsAndIDs.csv", sep = "/")
write.table(x =  input_files, file = filename, quote = FALSE, sep = ",",
            row.names = FALSE, col.names = FALSE)

