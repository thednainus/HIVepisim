#get values from ECDC model to estimate incidence rate

#parameters of my equation
#beta = transmission rate
#Incidence is proportional to beta x I(t)
# I(t) = number of infected individuals at time t
# converting beta to transmission probability (trans.prob) per act rate (act.rate)
# m = average number of partners per year per individual
# beta = m(1 -  (1 - trans.prob)^act.rate)

# after meeting with Erik on 27 January 2021
# We decided to fix act to 1.0 per day to estimate
# transmission probability per year
# we used mean partners = 1.29 (from paper Weiss et al 2020 (Epidemics))
# m= (mean_degree * 12 months)/duration = (1.29 * 12)/0.79
# m will be fixed at 19.59494 (based on equations above)
# act.rate in days * duration in days = 1 * 23.7 = 23.7 days
# 23.7 act rates per year per partner
act.rate <- 23.7
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

saveRDS(data_ecdc, "inst/data/transmission_probability_v2.RDS")
