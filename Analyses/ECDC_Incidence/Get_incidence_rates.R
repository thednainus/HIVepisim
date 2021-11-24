#get values from ECDC model to estimate incidence rate

#parameters of my equation
#beta = transmission rate
#Incidence is proportional to beta x I(t)
# I(t) = number of infected individuals at time t
# converting beta to transmission probability (trans.prob) per act rate (act.rate)
# m = average number of partners per year per individual
# beta = m(1 -  (1 - trans.prob)^act.rate)

# act.rate will be fixed at 81 acts per year (act rate for MSM reported for the
# last sexual partner in the USA: Wall et al 2013)
# 699 partners in past 3 months with median duration of 0.79 months.
#From Team meeting on 23 November 2021 with Erik
# we used mean partners = 1.29 (from paper Weiss et al 2020 (Epidemics))
# m= (mean_degree * 12 months)/duration = (1.29 * 12)/0.79
# m will be fixed at 17.77215 (based on equations above)
# act.rate in days * duration in days = 0.22 * 23.7 = 5.214 days
# act rate per year = 5.214
act.rate <- 5.214
m <- 19.59494

#N_Inf_M = annual number of newly-acquired HIV infections (ECDC model)
#N_Alive = total number of HIV-positive individuals who are still alive. This
#number is equal to the total number of infected individuals minus the observed
#number of individuals who dies minus the estimated number of individuals who
#died before being diagnosed with HIV

model1 <- readRDS("../../HIV_platform/Results_modelling1/HIVModelMainFit_20211022_110848.rds")


#create dataframe
data_ecdc <- data.frame(time = model1$Year,
                        incidence = model1$N_Inf_M,
                        infected_ind = model1$N_Alive)

# transmission probability per year
data_ecdc["trans.rate"] <- data_ecdc$incidence/data_ecdc$infected_ind

#1 - ((1 - (data_ecdc$trans.rate/m))^(1/act.rate))

data_ecdc["trans.prob.year"] <- 1 - ((1 - (data_ecdc$incidence)/(m * data_ecdc$infected_ind))^(1/act.rate))

# convert transmission probability per day
data_ecdc["trans.prob"] <- 1 - ((1 - data_ecdc$trans.prob.year)^(1/365))

saveRDS(data_ecdc, "inst/data/transmission_probability.RDS")
