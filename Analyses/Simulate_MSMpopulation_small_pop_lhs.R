# simulate MSM population based on lhs

library(EpiModel)
library(HIVepisim)
library(lubridate)
library(stringr)
library(lhs)

#get parameter values from file to run simulations

#for small pop size of npop1 = 1000 and npop2= npop1*3
#and inital pop size for pop2 was 300
#params1 <- sampler(n=1000, paramdf = data.frame(init_pop1_param = c(1, 300),
#                                                act.rate_param = c(0.1, 2),
#                                                inf.prob.param = c(0.005, 0.1)))

#for small pop size of npop1 = 5000 and npop2= npop1*5
#and initial pop size for pop2 was 500
#params3 <- sampler(n=1000, paramdf = data.frame(init_pop1_param = c(1, 2500),
#                                                act.rate_param = c(0.1, 2),
#                                                inf.prob.param = c(0.005, 0.1)))

#for small pop size of npop1 = 5000 and npop2= npop1*5
#and inital pop size for pop2 was 500
params4 <- sampler(n=1000, paramdf = data.frame(init_pop1_param = c(1, 650),
                                                act.rate_param = c(0.1, 2),
                                                inf.prob.param = c(0.005, 0.1)))

#params2 <- sampler(n=1000, paramdf = data.frame(init_pop1_param = c(50, 5000),
#                                                act.rate_param = c(0.1, 2),
#                                                inf.prob.param = c(0.005, 0.1)))

#saveRDS(params1, "inst/data/lhs_params_smallpop_20Sept2021.RDS")
#saveRDS(params2, "inst/data/lhs_params_largepop_20Sept2021.RDS")
#saveRDS(params3, "inst/data/lhs_params_smallpop_22Sept2021.RDS")
#saveRDS(params4, "inst/data/lhs_params_smallpop_23Sept2021.RDS")

line_number <- commandArgs(trailingOnly = TRUE)

params <- readRDS(system.file("data/lhs_params_smallpop_20Sept2021.RDS", package = "HIVepisim"))

params <- params[line_number,]

#npop1 <- 50000
#npop2 <- npop1 * 10

npop1 <- 5000
npop2 <- npop1 * 5

sim <- simulate_population(npop1 = npop1,
                           npop2 = npop2,
                           init_infect_pop2 = 500,
                           init.pop1.param = params$init_pop1_param,
                           act.rate.param = params$act.rate_param,
                           inf.prob.param = params$inf.prob.param)

saveRDS(sim, "results_sim.RDS")
