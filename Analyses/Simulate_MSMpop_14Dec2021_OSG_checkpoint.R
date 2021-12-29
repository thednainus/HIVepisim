# simulate MSM population based on lhs
# In here I only sampled the probability of initiating treatment once an HIV
# positive is diagnosed

# I also used the infection probability based on Erik's suggestion of getting
# this probability based on incidence in San Diego calculated using the ECD model

#param_list <- commandArgs(trailingOnly = TRUE)
#line_number <- as.numeric(param_list[1])
#message(line_number)

line_number <- 6607

#seed_value <- as.numeric(param_list[2])
#message(seed_value)
seed_value <- 60777

set.seed(seed_value)

#start of script
start_time <- Sys.time()

library(EpiModel)
library(EpiModelHPC)
library(HIVepisim)
library(HIVepisimAnalysis)
library(lubridate)
library(stringr)


#library(lhs)
#library(plyr)


# number of years to simulate
#years = 40

##beginning of simulation date
init_sim_date <- ymd("1980-01-01")
# end of simulation date
end_sim_date <- ymd("2021-01-01")

#number of days to simulate
total_steps <- as.numeric(end_sim_date - init_sim_date)

#total_steps <- 100

#checkpoint_steps
#checkpoint_steps <- 1872
checkpoint_steps <- total_steps
#checkpoint_steps <- 2

#if file data/sim1/sim1.cp.rda does not exist, it means that we are running this
# code for the first time
if(file.exists("data/sim1/sim1.cp.rda") == FALSE){

  #start counting
  iter <-  1

  saveRDS(iter, "iter.RDS")

  checkpoint_stage <- checkpoint_steps

  #just to check where I am at the checkpoint stage
  saveRDS(checkpoint_stage, "checkpoint_stage.RDS")


  # 63% per year is 0.0027 per day
  # 0.01253765 is 99% per year
  params4dim <- readRDS(system.file("data/params4dim_13Dec2021.RDS",
                                    package = "HIVepisim"))[line_number,]

  # Network Initialization --------------------------------------------------

  # time unit related to 1 day
  time.unit <- 1

  # Range of possible ages
  ages <- 18:80

  # Age-specific mortality rates for MALES for the San Diego
  # Values were obtained using https://wonder.cdc.gov/
  # group by County, Gender, Year, Age Group
  # and searching fro San Diego only
  # rates are per 100,000 population

  departure_rate_years <- c(1980:2016)

  #get mortality data for MALES for San Diego from 1979-1998
  mortality_1979_1998 <- read.delim(system.file("data/CompressedMortality_1979-1998.txt", package = "HIVepisim"))
  mortality_1999_2016 <- read.delim(system.file("data/CompressedMortality_1999-2016.txt", package = "HIVepisim"))

  dep1 <- list()
  years1 <- c(1980:1998)
  years2 <- c(1999:2016)

  for (i in 1:length(years1)){

    dep_year <- subset(mortality_1979_1998, Year == years1[i])$Crude.Rate[1:8]
    dep1[[i]] <- dep_year
  }

  for(j in 1:length(years2)){

    dep_year <- subset(mortality_1999_2016, Year == years2[j])$Crude.Rate[1:8]
    dep1[[j+19]] <- dep_year


  }

  names(dep1) <- c(paste("dep", years1, sep = ""), paste("dep", years2, sep = ""))

  departure_rates <- data.frame(dep1)


  # Per-capita daily death rate
  #dr_pp_pd <- lapply(departure_rates, function(x) x/100000)
  #dr_pp_pd <- lapply(dr_pp_pd, function(x) x/365)

  dr_pp_pd <- departure_rates/100000
  #assuming that 6.7% are MSM (Grey at al. 2016)
  dr_pp_pd <-  dr_pp_pd * 0.067
  dr_pp_pd <- dr_pp_pd/365


  # Create vector of daily death rates
  age_spans <- c(2, 5,  rep(10, 5), 6)
  #if using time unit per day
  #dr_vec <- rep(dr_pp_pd, times = age_spans)
  #if unsing time unit per week (7 days)
  dr_vec <- apply(dr_pp_pd, 2, function(x) rep(x, times = age_spans))
  dr_vec <- as.data.frame(dr_vec)
  #dr_vec["ages"] <- ages

  # Plot death rates
  #par(mar = c(3,3,2,1), mgp = c(2,1,0), mfrow = c(1,1))
  #plot(ages, dr_vec, type = "o", xlab = "age", ylab = "Mortality Rate")

  # Initialize network
  n_pop1 = 1000
  #n_pop1 = 30000
  n_pop2 = n_pop1 * 3

  #total number of individuals in the network
  n_total = n_pop1 + n_pop2
  nw <- network_initialize(n_total)

  # Set age attribute
  ageVec <- sample(ages, n_total, replace = TRUE)
  length(ageVec)
  nw <- set_vertex_attribute(nw, "age", ageVec)


  # Set origin attribute
  originVec1 <- rep("region", n_pop1)
  originVec2 <- rep("global", n_pop2)
  originVec <- c(originVec1, originVec2)
  print(table(originVec))


  nw <- set_vertex_attribute(nw, "origin", originVec)

  # Create vector of diagnose statuses
  diagStatusVec1 <- rep(0, n_pop1)

  #n_inf_pop1 <- 300.9307
  n_inf_pop1 <- params4dim$init_pop1_param
  init.Infected1 <- sample(1:n_pop1, n_inf_pop1)
  diagStatusVec1[init.Infected1] <- 1


  diagStatusVec2 <- rep(0, n_pop2)
  #n_inf_pop2 <- 250
  n_inf_pop2 <- n_pop2 * 0.03
  init.Infected2 <- sample(1:n_pop2, n_inf_pop2)
  diagStatusVec2[init.Infected2] <- 1

  diagStatusVec <- c(diagStatusVec1, diagStatusVec2)
  print(table(diagStatusVec))

  # Set status attribute
  nw <- set_vertex_attribute(nw, "diag.status", diagStatusVec)
  print(nw)

  # set vector of status
  statusVector <- diagStatusVec
  statusVector[diagStatusVec == 1] <- "i"
  statusVector[diagStatusVec == 0] <- "s"

  print(table(statusVector))

  # Set status attribute
  nw <- set_vertex_attribute(nw, "status", statusVector)
  #nw

  # Create vector of migrant status
  migrantVec1 <- rep(1, n_pop1)
  migrantVec2 <- rep(2, n_pop2)

  migrantVec <- c(migrantVec1, migrantVec2)
  print(table(migrantVec))

  # Set status attribute
  nw <- set_vertex_attribute(nw, "migrant", migrantVec)
  #nw



  # Create vector of risk groups
  # 20% of population will be on risk group 2
  riskGroupVector1 <- sample(x = 1:2, size = n_pop1, replace = TRUE, prob = c(0.8, 0.2))
  riskGroupVector2 <- sample(x = 1:2, size = n_pop2, replace = TRUE, prob = c(0.8, 0.2))

  #riskGroupVector1 <- rep(1:2, each = n_pop1/2)
  #riskGroupVector2 <- rep(1:2, each = n_pop2/2)

  riskGroupVector <- c(riskGroupVector1, riskGroupVector2)
  print(table(riskGroupVector))

  # Set risk group attribute
  nw <- set_vertex_attribute(nw, "risk.group", riskGroupVector)
  print(nw)




  # Network Model Estimation ------------------------------------------------

  # Review network object
  print(nw)

  # Define the formation model: edges

  # Density (edges) sets the baseline probability of forming a tie, and it is
  # tipically included in all models.

  # Heterogeneity by attribute (nodefactor) allows the probability of an edge to
  # depend on the nodal attribute "risk.group", for example.

  # Selective mixing by attribute (nodematch) allows the probability of an edge
  # to depend on the attributes of both nodes. It is used to capture the level of
  # assortative (or disassortative) mixing between groups.

  # Degree distribution (concurrent): There are many ways to specify further
  # heterogeneity in the degree distribution. In the absence of further
  # specification, the conditional probability of a partnership forming is
  # Bernoulli with the (group-specific) parameter determined by coefficients
  # on the terms above, and the resulting degree distribution is a binomial
  # mixture. Here we will modify this by specifying the number of nodes
  # which have more than one partnership at any time.


  #formation <- ~edges + nodemix("origin") + nodefactor("risk.group") + nodematch("risk.group") + concurrent
  # nodemix as specified below will only estimate edges for global.global and region.region
  # to know the other of terms you can type
  # summary(nw ~ nodemix("origin"))
  formation <- ~edges + nodemix("origin", levels2 = c(1, 2))

  #formation <- ~edges + nodemix("origin", levels2 = c(1, 2)) + concurrent(by="origin")


  # CHECK how to specify degree and concurrent for both populations
  #formation <- ~edges + nodemix("origin", levels2 = c(1, 2)) + degree(0, by = "origin") +
  #  concurrent(by="origin")
  #formation <- ~edges

  # Target statistics
  #target.stats <- c(250, 375, 225, 100)
  # target.stats below is including nodemix for global.global and region.region
  #overall edges = 250; edges b/t global.global = 200; edges b/t global.region = 0
  # then edges b/t global.region.region = 250 - 150 - 0 = 100
  #target.stats <- c((0.04 * n_pop1/2 + 0.04 * n_pop2/2), 0.04 * n_pop2/2, 0)
  # mean degree of 1.17 from Weiss et al 2020 (Epidemics)
  target.stats <- c((1.19 * n_pop1/2 + 1.19 * n_pop2/2), 1.19 * n_pop2/2, 0)
  #target.stats <- c((1.19 * n_pop1/2 + 1.19 * n_pop2/2), 1.19 * n_pop2/2, 0,
  #                  0.263 * n_pop2, 0.263 * n_pop1)

  #target.stats <- c((0.2 * n_pop1/2 + 0.2 * n_pop2/2),
  #                  0.2 * n_pop2/2,
  #                  0,
  #                  0.8 * n_pop2,
  #                  0.8 * n_pop1,
  #                  0,
  #                  0)

  print(target.stats)

  # Dissolution  model
  # 200 is the time step for a dissolution to happen
  #coef.diss <- dissolution_coefs(~offset(edges), 200, mean(dr_vec))
  #durs <- c(60, 60, 2)
  durs <- c(23.7, 23.7, 2)
  coef.diss <- dissolution_coefs(~offset(edges) +
                                 offset(nodemix("origin", levels2 = c(1, 2))),
                                 duration = durs, d.rate = mean(unlist(lapply(dr_vec, mean))))
  print(coef.diss)

  # Fit the model
  # # Fit the TERGM
  # to simulate large networks

  if(n_pop1 > 10000){

    est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE,
                  set.control.ergm = control.ergm(MCMC.burnin=1e7, MCMC.interval=1e7,
                                                  MCMC.samplesize=20000,
                                                  init.MPLE.samplesize = 1e8,
                                                  SAN = control.san(SAN.nsteps = 1e8)
                  ))

    print(est)

    # to simulate large networks
    dx <- netdx(est, nsims = 5, nsteps = 500, keep.tedgelist = TRUE,
                set.control.stergm = control.simulate.network(MCMC.maxchanges = 1e7,
                                                              MCMC.burnin.min= 1.5e5,
                                                              MCMC.burnin.max =1.5e5))

  } else {
    #est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE)
    est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE)
    dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE,
                set.control.stergm = control.simulate.network(MCMC.maxchanges = 1e7)
                )

    # Model diagnostics
    # Simulate time series to examine timed edgelist

    # You can ignore message:
    #   Warning message:
    #   In x$coef.diss$duration^2 * dgeom(2:(nsteps + 1), 1/x$coef.diss$duration) :
    #   longer object length is not a multiple of shorter object length
    #
    # https://github.com/statnet/EpiModel/issues/529
  }


  #dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE,
  #            nwstats.formula = ~edges + nodemix("origin", levels2 = c(1, 2)) + degree(0:3, by = "origin"))
  print(est)
  print(dx)
  save(est, file = "est.rda")
  #plot(dx)
  # Extract timed-edgelist
  #te <- as.data.frame(dx)
  #head(te)

  # Limit to non-censored edges
  #te2 <- te[which(te$onset.censored == FALSE & te$terminus.censored == FALSE),
  #          c("head", "tail", "duration")]
  #head(te2)

  # Look up the age group of head and tail nodes
  #te2$ag.head <- originVec[te2$head]
  #te2$ag.tail <- originVec[te2$tail]
  #head(te2)

  # Calculate mean durations in each age-mixed group
  #mean(te$duration[te2$ag.head == "global" & te2$ag.tail == "global"])
  #mean(te$duration[te2$ag.head != te2$ag.tail])
  #mean(te$duration[te2$ag.head == "region" & te2$ag.tail == "region"])

  # Compare to targets
  #durs



  # EpiModel Model Simulation -----------------------------------------------
  # Base model
  #a1.rate = 3.798435e-05,
  #a2.rate = 3.798435e-05 * 0.0645
  #m21.rate = 0.00129

  # parameters from epimodelHIV param_msm
  # hiv.test.rate = 0.01325 per week. Per day will be 0.01325/7
  # test.window.int = 21/7 per week. Per day will be 21
  # tx.init.prob = 0.092 per week, tx.init.prob = 0.092/7 per day
  # tx.halt.prob = 0.0102 per week, tx.halt.prob = 0.0102/7 per day
  # tx.reinit.prob = 0.00066 per week, tx.reinit.prob = 0.00066 per day
  # a1.rate = 0.00052 per week, a1.rate = 0.00052/7 per day
  # a2.rate = 0.00052 per week, a2.rate = 0.00052/7 per day
  # a2.rate = 0.00005 per week
  # art_start = 24 * 365

  #in weeks
  #art_start <- 24 * 52
  #in days
  #in case there is checkpoint, so we can read the paramenters values again
  save(params4dim, time.unit, departure_rate_years, dr_vec, file = "params.rda")
}

if(file.exists("data/sim1/sim1.cp.rda") == TRUE){

  #in case checkpoint exists
  #params2dim <- readRDS("params2dim.RDS")
  load("params.rda")
  iter = readRDS("iter.RDS")
  checkpoint_stage <- readRDS("checkpoint_stage.RDS")

}


year_art_start <- ymd("2004-01-01")
art_start <- as.numeric(year_art_start - init_sim_date)


#hiv.test.rate (mean probability of HIV testing per day for MSM)
# parameter was fixed to 0.045
# see notes on my notebook (pages 95-96 for reason)
# probably not the right reason. Think a better way of getting this parameter
# Age-specific mortality rates for MALES for the San Diego
# Values were obtained using https://wonder.cdc.gov/
# group by County, Gender, Year, Age Group
# and searching fro San Diego only
# rates are per 100,000 population
a1.rate.times <- c(2000:2019)
#crude birth rate based on crude birth rate table for San Diego County
#this is all births (males + females) from 2000 to 2019

a1.rates <- c(44272/2813833, 43758/2849238, 43951/2890256,
              45368/2927216, 45758/2953703, 45897/2966783,
              46876/2976492, 47545/2998477, 46742/3032689,
              44960/3064436, 44838/3095314, 43621/3125265,
              44391/3161751, 43627/3201418, 44596/3235143,
              43960/3267933, 42741/3287280, 40889/3309627,
              39921/3333127, 38445/3351784)
#assuming that 50% are males
#a1.rates <- a1.rates * 0.5
#assuming that 6.7% are MSM (Grey at al. 2016)
a1.rates <- a1.rates * 0.067
#now convert rate per year to rate per day
a1.rates <- a1.rates/365

#increase arrival rate because by the end of 40 years not many individuals
#arrive at the network
a1.rates <- a1.rates * 10




a2.rate.times <- c(1980:2019)
#crude birth rate based on
#https://data.worldbank.org/indicator/SP.DYN.CBRT.IN?end=2019&start=1980&view=chart

a2.rates <- c(27.421, 27.892, 28.092, 27.487, 27.222, 27.254, 27.332,
              27.275, 26.79, 26.274, 25.879, 25.214, 24.582, 24.182,
              23.804, 23.368, 23.079, 22.75, 22.337, 21.912, 21.677,
              21.329, 21.071, 20.859, 20.703, 20.57, 20.422, 20.341,
              20.232, 20.036, 19.809, 19.629, 19.514, 19.304, 19.217,
              18.965, 18.949, 18.629, 18.169, 17.897)
#rates are reported to 1,000 individuals
a2.rates <- a2.rates/1000
#assuming that 50% are males
#a2.rates <- a2.rates * 0.5
#assuming that 6.7% are also MSM
a2.rates <- a2.rates * 0.067
#now convert rate per year to rate per day
a2.rates <- a2.rates/365

#increase arrival rate because by the end of 40 years not many individuals
#arrive at the network
a2.rates <- a2.rates * 2

#migration in both directions will be set to 1 individual in 10 years

#hiv.test.rate:data was previously fit using a fixed value of 0.045
hiv.test.rate.df <- readRDS(system.file("data/probability_msm_tested_per_day.RDS",
                                        package = "HIVepisim"))

#tx.init.prob = 0.0027 (check reason on my notes page 136)

inf.prob.param <-  readRDS(system.file("data/transmission_probability.RDS",
                                       package = "HIVepisim"))
inf.prob.param$trans.prob <- inf.prob.param$trans.prob * params4dim$multiplier
#act.rate.param <- 0.22
act.rate.param <- params4dim$act.rate

# 63% per year is 0.0027 per day
# 0.01253765 is 99% per year
#tx.init.prob.param <- runif(10000, min = 0.0001, max = 0.01253765)
#saveRDS(tx.init.prob.param, "inst/data/tx.init.prob.param.RDS")
#line_number <- 1
#get line number to simulate MSM population
#line_number <- commandArgs(trailingOnly = TRUE)
#line_number <- as.numeric(line_number)
#print(line_number)
#tx.init.prob.param <- readRDS(system.file("data/tx.init.prob.param.RDS",
#                                          package = "HIVepisim"))[line_number]

#rates to progress to each stage is the average time for all individuals to
#transition from one CD4+ category to the next (Cori et al. 2015)
# analyses were reported for 1,039 individuals
#rate per individual were then divided by 1039 (Cori et al. 2015)
param <- param.net(time.unit = time.unit,
                   init_date = init_sim_date,
                   groups = 1,
                   act.rate = act.rate.param,
                   stage_prog_rate0 = (1/((0.5 * 365) / time.unit))/1039,
                   stage_prog_rate1 = (1/((3.32 * 365) / time.unit))/1039,
                   stage_prog_rate2 = (1/((2.7 * 365) / time.unit))/1039,
                   stage_prog_rate3 = (1/((5.5 * 365) / time.unit))/1039,
                   f1 = 0.76,
                   f2 = 0.19,
                   f3 = 0.05,
                   f4 = 0,
                   hiv.test.rate = hiv.test.rate.df,
                   test.window.int = 21 / time.unit,
                   art_start = art_start,
                   tx.init.prob = params4dim$tx.init.prob.param,
                   tx.halt.prob = 0,
                   tx.reinit.prob = 0,
                   inf.prob = inf.prob.param,
                   ws0 = 1,
                   ws1 = 0.1,
                   ws2 = 0.1,
                   ws3 = 0.1,
                   ws4 = 0.3,
                   wc1 = 1,
                   wc2 = 0.5,
                   wc3 = 0.05,
                   wr1 = 1,
                   wr2 = 10,
                   aids.mr = (1/((5.06 * 365) / time.unit))/1039,
                   asmr = list(dr_times = departure_rate_years, dep_vec = dr_vec),
                   a1.rate = list(a1.times = a1.rate.times, a1.rates = a1.rates),
                   a2.rate = list(a2.times = a2.rate.times, a2.rates = a2.rates),
                   arrival.age = 18,
                   m12.rate = 1/((n_pop1/50)*365),
                   m21.rate = 1/((n_pop2/50)*365))
# In order to have a m12.rate that will generate 500 migrations every 10 years,
# the m12.rate should be 1/(100 * 365) (rate per day)
# to balance number of migrant for m21.rate, we will have 15000/x = 50 indvidulas per year;
# x = 15000/50; x = 300 years
# m21.rate = 1/(300*365)


init <- init.net()

#in weeks
#nsteps = years * 52
# in days

control <- control.net(type = NULL, simno = 1, currsim = 1,
                       nsteps = checkpoint_stage, start = 1, nsims = 1,
                       when2save_at = total_steps,
                       ncores = 1, resimulate.network = TRUE, tergmLite = TRUE,
                       initialize.FUN = initialize_mig,
                       resim_nets.FUN = resim_nets,
                       hivtest.FUN = hivtest_msm,
                       hivtx.FUN = hivtx_msm,
                       hivprogress.FUN = hivprogress_msm,
                       hivtrans.FUN = hivtrans_mig,
                       aging.FUN = aging_msm,
                       departure.FUN = departure_mig,
                       migration.FUN = migration,
                       arrivals.FUN = arrivals_mig,
                       nwupdate.FUN = nwupdate_mig,
                       prevalence.FUN = prevalence_mig,
                       verbose.FUN = verbose.net,
                       save.nwstats = TRUE,
                       save_nodes = TRUE,
                       nwstats.formula = ~edges + nodemix("origin", levels2 = c(1, 2)) + nodefactor("origin", levels = c("global", "region")) + concurrent,
                       save.transmat = TRUE,
                       save.stats = TRUE,
                       verbose = TRUE,
                       )


if(iter < 64){



  # I added a + 1 in  nsteps = total_steps + 1
  #because I don't want to save the other files I generated
  # at the very end of the simulations

  #I added the when2save_at, so I can introduce checkpoint for OSG
  #terminate files when I reach checkpoint_steps
  #but only save all files at the very end



# creating stop condition to exit with code 85
  #sim <- netsim_hpc("est.rda", param, init, control,
  #                  cp.save.int = checkpoint_steps, save.min=FALSE, save.max=FALSE,
  #                  compress = "gzip")
  start_time <- Sys.time()
  sim <- netsim(est, param, init, control)

  #end of script
  end_time <- Sys.time()
  print("Simulation took:")
  print(end_time - start_time)

  iter = iter + 1
  saveRDS(iter, "iter.RDS")

  checkpoint_stage <- checkpoint_stage + checkpoint_steps
  saveRDS(checkpoint_stage, "checkpoint_stage.RDS")

  #quit with exit code 85
  #so condor can create the checkpoint in SPOOL
  #print(quit(save = "no", status = 85))
  #quit(save = "no", status = 85)
}


if(iter == 64){

  sim <- netsim_hpc("est.rda", param, init, control,
                    cp.save.int = checkpoint_steps, save.min=FALSE, save.max=FALSE,
                    compress = "gzip")
  #quit(status = 0)
}







saveRDS(sim, "results_sim.RDS")

#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

sim_df <- as.data.frame(sim)
#sim_df$sim <- "1"
sim_df$type <- paste("param", as.character(line_number), sep = "_")

pop1_data <- sim_df[c("type", "time", "newDx_pop1", "incid.pop1", "num.pop1", "i.num.pop1")]
pop1_data_freq <- data.frame(type = pop1_data["type"],
                             time = pop1_data["time"],
                             newDx_pop1 = pop1_data["newDx_pop1"],
                             incidence_pop1 = pop1_data["incid.pop1"])


pop1_data_freq["time_years"] <- days2years(pop1_data_freq$time, init_date = init_sim_date)

pop1_data <- pop1_data_freq
pop1_data["year"] <- unlist(lapply(pop1_data$time_years, function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
pop1_data$year <- as.factor(pop1_data$year)
pop1_data$type <- as.factor(pop1_data$type)
#pop1_data$sim <- as.factor(pop1_data$sim)
pop1_data <- pop1_data[c(6, 1:4)]

#aggregate by new HIV diagnosis and by incidence for population 1
newDx_pop1_agg <- aggregate(newDx_pop1  ~ year + type, data = pop1_data, FUN=sum)
incid_pop1_agg <- aggregate(incid.pop1  ~ year + type, data = pop1_data, FUN=sum)


newDx_pop1_agg["rep_param"] <- newDx_pop1_agg$type
incid_pop1_agg["rep_param"] <- incid_pop1_agg$type

#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
incidence_model <- readRDS(system.file("data/ECDC_incidence_model_22Oct2021.RDS",
                                       package = "HIVepisimAnalysis"))



#new HIV diagnosis
log_weights_dx <- compute_log_importance_weight_newDx(incidenceDiag$frequency,
                                                      diag_sim = newDx_pop1_agg)
log_weights_dx_df <- data.frame(rep_param = newDx_pop1_agg$rep_param[1],
                                log_weights = log_weights_dx)

#new HIV incidence
log_weights_incid <- compute_log_importance_weight_incidence(incidence_model$ECDC_incidence.N_Inf_M, incid_sim = incid_pop1_agg[1:41,])
log_weights_incid_df <- data.frame(rep_param = incid_pop1_agg$rep_param[1], log_weights = log_weights_incid)


saveRDS(log_weights_dx_df, "log_weights_dx_df.RDS")
saveRDS(newDx_pop1_agg, "summary_newDx_pop1.RDS")


saveRDS(log_weights_incid_df, "log_weights_incid_df.RDS")
saveRDS(incid_pop1_agg, "summary_incidence_pop1.RDS")
quit(save = "no", status = 0)

