library(EpiModel)
library(HIVepisim)
library(EpiModelHPC)

rm(list = ls())

# PARAMETERS are tests and should be CONFIRMED

# Network Initialization --------------------------------------------------

# time unit related to 1 day
time.unit <- 1

# Range of possible ages
ages <- 18:80

# numner of years to simulate
years = 5


# Age-specific mortality rates for MALES for the UK in 2018
# (https://www.statista.com/statistics/1125118/death-rate-united-kingdom-uk-by-age/)
# Rates per 1,000 for age groups: 15-19, 20-24, 25-29,
#                                   30-34, 35-39, 40-44, 45-49, 50-54, 55-59,
#                                   60-64, 65-69, 70-74, 75-79, 80-84
departure_rate <- c(0.3, 0.5, 0.6, 0.8, 1.2, 1.8, 2.7, 3.9, 5.8, 9.1, 14.6, 22.5,
                    39.2, 68.2)

# This is just to test the code so departures can happen faster
#departure_rate <- departure_rate * 0.2

# Per-capita daily death rate
dr_pp_pd <- departure_rate / 1000 / 365
dr_pp_pw <- dr_pp_pd * time.unit


# Create vector of daily death rates
age_spans <- c(2, rep(5, 12), 1)
#if using time unit per day
#dr_vec <- rep(dr_pp_pd, times = age_spans)
#if unsing time unit per week (7 days)
dr_vec <- rep(dr_pp_pw, times = age_spans)
data.frame(ages, dr_vec)

# Plot death rates
#par(mar = c(3,3,2,1), mgp = c(2,1,0), mfrow = c(1,1))
#plot(ages, dr_vec, type = "o", xlab = "age", ylab = "Mortality Rate")

# Initialize network
n_pop1 = 500
n_pop2 = 500
#n_pop1 = 10
#n_pop2 = 10
n_total = n_pop1 + n_pop2
#n <- 50000
nw <- network_initialize(n_total)

# Set age attribute
ageVec <- sample(ages, n_total, replace = TRUE)
length(ageVec)
nw <- set_vertex_attribute(nw, "age", ageVec)


# Set origin attribute
originVec1 <- rep("region", n_pop1)
originVec2 <- rep("global", n_pop2)
originVec <- c(originVec1, originVec2)
table(originVec)


nw <- set_vertex_attribute(nw, "origin", originVec)

# Create vector of diagnose statuses
diagStatusVec1 <- rep(0, n_pop1)
init.Infected1 <- sample(1:n_pop1, 0.005 * n_pop1)
#init.Infected1 <- sample(1:n_pop1, 0.4 * n_pop1)
diagStatusVec1[init.Infected1] <- 1


diagStatusVec2 <- rep(0, n_pop2)
init.Infected2 <- sample(1:n_pop2, 0.015 * n_pop2)
#init.Infected2 <- sample(1:n_pop2, 0.4 * n_pop2)
diagStatusVec2[init.Infected2] <- 1

diagStatusVec <- c(diagStatusVec1, diagStatusVec2)
table(diagStatusVec)

# Set status attribute
nw <- set_vertex_attribute(nw, "diag.status", diagStatusVec)
nw

# set vector of status
statusVector <- diagStatusVec
statusVector[diagStatusVec == 1] <- "i"
statusVector[diagStatusVec == 0] <- "s"

table(statusVector)

# Set status attribute
nw <- set_vertex_attribute(nw, "status", statusVector)
nw

# Create vector of migrant status
migrantVec1 <- rep(1, n_pop1)
migrantVec2 <- rep(2, n_pop2)

migrantVec <- c(migrantVec1, migrantVec2)
table(migrantVec)

# Set status attribute
nw <- set_vertex_attribute(nw, "migrant", migrantVec)
nw



# Create vector of risk groups
riskGroupVector1 <- rep(1:2, each = n_pop1/2)
riskGroupVector2 <- rep(1:2, each = n_pop1/2)

riskGroupVector <- c(riskGroupVector1, riskGroupVector2)
table(riskGroupVector)

# Set risk group attribute
nw <- set_vertex_attribute(nw, "risk.group", riskGroupVector)
nw




# Network Model Estimation ------------------------------------------------

# Review network object
nw

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
#formation <- ~edges

# Target statistics
#target.stats <- c(250, 375, 225, 100)
# target.stats below is including nodemix for global.global and region.region
#overall edges = 250; edges b/t global.global = 200; edges b/t global.region = 0
# then edges b/t global.region.region = 250 - 150 - 0 = 100
target.stats <- c((0.04 * n_pop1/2 + 0.04 * n_pop2/2), 0.04 * n_pop1/2, 0)
#target.stats <- c((0.04 * n_pop1/2 + 0.04 * n_pop2/2),
#                  0.04 * n_pop1/2,
#                  0,
#                  0.538 * n_pop1,
#                  0.538 * n_pop2)
#target.stats <- c(9, 4, 0)
# target.stat below is for ~edges only
target.stats

# Dissolution  model
# 200 is the time step for a dissolution to happen
#coef.diss <- dissolution_coefs(~offset(edges), 200, mean(dr_vec))
#durs <- c(60, 60, 2)
durs <- c(23.7, 23.7, 2)
coef.diss <- dissolution_coefs(~offset(edges) +
                                 offset(nodemix("origin", levels2 = c(1, 2))),
                               duration = durs, d.rate = mean(dr_vec))
coef.diss

# Fit the model
# # Fit the TERGM
est <- netest(nw, formation, target.stats, coef.diss)
save(est, file = "fit.rda")

# Model diagnostics
# Simulate time series to examine timed edgelist
dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE)

# Extract timed-edgelist
te <- as.data.frame(dx)
head(te)

# Limit to non-censored edges
te2 <- te[which(te$onset.censored == FALSE & te$terminus.censored == FALSE),
          c("head", "tail", "duration")]
head(te2)

# Look up the age group of head and tail nodes
te2$ag.head <- originVec[te2$head]
te2$ag.tail <- originVec[te2$tail]
head(te2)

# Calculate mean durations in each age-mixed group
mean(te$duration[te2$ag.head == "global" & te2$ag.tail == "global"])
mean(te$duration[te2$ag.head != te2$ag.tail])
mean(te$duration[te2$ag.head == "region" & te2$ag.tail == "region"])

# Compare to targets
durs



# EpiModel Model Simulation -----------------------------------------------
# Base model
#a1.rate = 3.798435e-05,
#a2.rate = 3.798435e-05 * 0.0645
#m21.rate = 0.00129
param <- param.net(time.unit = time.unit,
                   groups = 1,
                   stage_prog_rate1 = 1/((3.32 * 365) / time.unit),
                   stage_prog_rate2 = 1/((2.7 * 365) / time.unit),
                   stage_prog_rate3 = 1/((5.5 * 365) / time.unit),
                   f1 = 0.76,
                   f2 = 0.19,
                   f3 = 0.05,
                   f4 = 0,
                   hiv.test.rate = 0.01325,
                   test.window.int = 21/7,
                   tx.init.prob = 0.092,
                   tx.halt.prob = 0.0102,
                   tx.reinit.prob = 0.00066,
                   trans.r =5e-05,
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
                   aids.mr = 1/((5.06 * 365) / time.unit),
                   asmr = dr_vec,
                   a1.rate = 0.00052,
                   a2.rate = 0.00052,
                   arrival.age = 18,
                   m12.rate = 0.00016125,
                   m21.rate = 0.00016125)

init <- init.net()
#6752
control <- control.net(type = NULL, nsteps = 1 + (365 * years), start = 1,nsims = 1,
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
                       save.nwstats = FALSE,
                       save.transmat = TRUE,
                       verbose = TRUE)

sim <- netsim(est, param, init, control)

plot(sim, y = c("i.num.pop1", "i.num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = "i.num.pop2", qnts = 1, legend = TRUE)
plot(sim, y = "i.num.pop1", qnts = 1, legend = TRUE)
plot(sim, y = c("num.pop1", "num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("s.num.pop1", "s.num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = "s.num.pop2", qnts = 1, legend = TRUE)
plot(sim, y = "s.num.pop1", qnts = 1, legend = TRUE)
plot(sim, y = c("i.prev.pop1", "i.prev.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = "i.prev.pop2", qnts = 1, legend = TRUE)
plot(sim, y = "i.prev.pop1", qnts = 1, legend = TRUE)
plot(sim, y = "nArrivals_mig1", qnts = 1, legend = TRUE)
plot(sim, y = "nArrivals_mig2", qnts = 1, legend = TRUE)




quartz()

par(mfrow = c(1, 1))
plot(sim, y = "i.num.pop2", qnts = 1, main = "Total Prevalence", mean.col = 3, qnts.col = 3)
plot(sim_notergmLite, y = "i.num.pop2", qnts = 1, mean.col = 6, qnts.col = 6, add = TRUE)
legend("topleft", c("tergmLite", "No tergmLite"), lwd = 3, col = c(3, 6), bty = "n")



print(date())
sim <- netsim(est, param, init, control)
print(date())

sim <- netsim_hpc(x = "fit.rda", param = param, init = init, control = control,
                   cp.save.int = 100, save.min = TRUE,save.max = TRUE,
                   compress = TRUE, verbose = TRUE)

sim <- netsim_hpc(x = "data/sim1/sim.cp.rda", param = param, init = init, control = control,
                  cp.save.int = 100, save.min = TRUE,save.max = TRUE,
                  compress = TRUE, verbose = TRUE)

load("data/sim1/sim.cp.rda")


# Simulation plot ----
plot(sim, qnts = 1)

# Mean age summary statistic
plot(sim, y = "age.mean")

# Export data to data frame
df <- as.data.frame(sim, out = "mean")
head(df$age.mean)
tail(df$age.mean)

# Population size over time
plot(sim, y = "num")

# Deaths per day
plot(sim, y = c("hstage1", "hstage2", "hstage3"), qnts = FALSE, legend = TRUE)
plot(sim, y = "hstage0")

plot(sim, y = "i.prev")
plot(sim, y = "nNew")
plot(sim, y = "incid")


# Transmission matrix ----
tm <- get_transmat(sim)
transphylo <- as.phylo.transmat(tm)
tmbytree <- get.transmat.phylo(tm)
# to make sure that transmission only happens between individuals from same origin
test <- tm[tm$infOrigin != tm$susOrigin,]

origin_inf <- read.csv("infected_origin.csv")

transphylo[[1]]$tip.label %in% origin_inf$infID

