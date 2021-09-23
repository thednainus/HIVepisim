# simulate MSM population based on lhs

library(EpiModel)
library(HIVepisim)
library(lubridate)
library(stringr)
library(lhs)

#get parameter values from file to run simulations
#params_list <- sampler(n=1000)
#saveRDS(params_list, "inst/data/lhs_params_17Sept2021.RDS")

#params1 <- sampler(n=1000, paramdf = data.frame(init_pop1_param = c(1, 300),
#                                                act.rate_param = c(0.1, 2),
#                                                inf.prob.param = c(0.005, 0.1)))

#params2 <- sampler(n=1000, paramdf = data.frame(init_pop1_param = c(50, 5000),
#                                                act.rate_param = c(0.1, 2),
#                                                inf.prob.param = c(0.005, 0.1)))

#saveRDS(params1, "inst/data/lhs_params_smallpop_20Sept2021.RDS")
#saveRDS(params2, "inst/data/lhs_params_largepop_20Sept2021.RDS")

line_number <- commandArgs(trailingOnly = TRUE)

params <- readRDS(system.file("data/lhs_params_largepop_20Sept2021.RDS", package = "HIVepisim"))

params <- params[line_number,]

npop1 <- 50000
npop2 <- npop1 * 10


sim <- simulate_population(npop1 = npop1,
                           npop2 = npop2,
                           init_infect_pop2 = 1000,
                           init.pop1.param = params$init_pop1_param,
                           act.rate.param = params$act.rate_param,
                           inf.prob.param = params$inf.prob.param)

saveRDS(sim, "results_sim.RDS")
