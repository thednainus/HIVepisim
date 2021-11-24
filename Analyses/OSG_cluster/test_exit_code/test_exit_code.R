#load("Analyses/OSG_cluster/test_exit_code/i_iter_test.rdata")
load("i_iter_test.rdata")
if(iter < 3){
  iter = iter + 1
  save(iter, file = "i_iter_test.rdata")
  quit(save = "no", status = 85)
}

if(iter == 3){
  print("I'm in the final if")
  quit(save = "no", status = 0)
}


