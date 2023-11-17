reps <- seq(1, 30, 1)
reps <- paste("rep", reps, sep = "_")

dir_names <- list("best_trajectories_250migrants",
                  "best_trajectories_500migrants",
                  "best_trajectories_750migrants")
params <- list("params_1067",
               "params_2348")

filename <- "processing_network_results.tar.gz"

for (dir_name in dir_names){

  for (param in params){

    for (i in 1:length(reps)){

      dir1 <- paste("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper",
                    dir_name, param, reps[i], filename, sep = "/")

      dir2 <- paste("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/deepseq",
                    dir_name, param, reps[i], sep = "/")

      if(dir.exists(dir2) == TRUE){
        commands <- paste("mv", dir1, dir2, sep = " ")
        print(commands)
        system(commands)
      }


    }

  }


}



