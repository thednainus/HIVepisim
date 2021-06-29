library(EpiModel)
library(HIVepisim)


rm(list = ls())

# PARAMETERS are tests and should be CONFIRMED

# Network Initialization --------------------------------------------------

# time unit related to 1 day
time.unit <- 1

# Range of possible ages
ages <- 18:80

# numner of years to simulate
years = 40


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
n_pop1 = 1000
n_pop2 = 2000
#n_pop1 = 50000
#n_pop2 = 950000
#n_pop1 = 5
#n_pop2 = 5
#n_pop1 = 50
#n_pop2 = 50
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
#n_inf_pop1 <- 50
n_inf_pop1 <- 10
init.Infected1 <- sample(1:n_pop1, n_inf_pop1)
diagStatusVec1[init.Infected1] <- 1


diagStatusVec2 <- rep(0, n_pop2)
#n_inf_pop2 <- 250
n_inf_pop2 <- 15
init.Infected2 <- sample(1:n_pop2, n_inf_pop2)
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
# 20% of population will be on risk group 2
riskGroupVector1 <- sample(x = 1:2, size = n_pop1, replace = TRUE, prob = c(0.8, 0.2))
riskGroupVector2 <- sample(x = 1:2, size = n_pop2, replace = TRUE, prob = c(0.8, 0.2))

#riskGroupVector1 <- rep(1:2, each = n_pop1/2)
#riskGroupVector2 <- rep(1:2, each = n_pop2/2)

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

# CHECK how to specify degree and concurrent for both populations
#formation <- ~edges + nodemix("origin", levels2 = c(1, 2)) + degree(0, by = "origin") +
#  concurrent(by="origin")
#formation <- ~edges

# Target statistics
#target.stats <- c(250, 375, 225, 100)
# target.stats below is including nodemix for global.global and region.region
#overall edges = 250; edges b/t global.global = 200; edges b/t global.region = 0
# then edges b/t global.region.region = 250 - 150 - 0 = 100
target.stats <- c((0.04 * n_pop1/2 + 0.04 * n_pop2/2), 0.04 * n_pop2/2, 0)

#target.stats <- c((0.2 * n_pop1/2 + 0.2 * n_pop2/2),
#                  0.2 * n_pop2/2,
#                  0,
#                  0.8 * n_pop2,
#                  0.8 * n_pop1,
#                  0,
#                  0)

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
# to simulate large networks

# est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE,
#               set.control.ergm = control.ergm(MCMC.burnin=1e8, MCMC.interval=1e8,
#                                               MCMC.samplesize=30000,
#                                               init.MPLE.samplesize = 1e9,
#                                               SAN.control = control.san(SAN.nsteps = 1e9)
#                                               ))


est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE)


#save(est, file = "fit.rda")

# Model diagnostics
# Simulate time series to examine timed edgelist

# You can ignore message:
#   Warning message:
#   In x$coef.diss$duration^2 * dgeom(2:(nsteps + 1), 1/x$coef.diss$duration) :
#   longer object length is not a multiple of shorter object length
#
# https://github.com/statnet/EpiModel/issues/529


# to simulate large networks
# dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE,
#             set.control.ergm = control.simulate.ergm(MCMC.init.maxchanges = 1e8,
#                                                      MCMC.burnin.min= 1.5e5,
#                                                      MCMC.burnin.max =1.5e5))



dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE)

#dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE,
#            nwstats.formula = ~edges + nodemix("origin", levels2 = c(1, 2)) + degree(0:3, by = "origin"))
#dx
plot(dx)
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
art_start <- 24 * 365
#art_start <- 14 * 365


param <- param.net(time.unit = time.unit,
                   groups = 1,
                   act.rate = 0.60,
                   stage_prog_rate0 = 1/((0.5 * 365) / time.unit),
                   stage_prog_rate1 = 1/((3.32 * 365) / time.unit),
                   stage_prog_rate2 = 1/((2.7 * 365) / time.unit),
                   stage_prog_rate3 = 1/((5.5 * 365) / time.unit),
                   f1 = 0.76,
                   f2 = 0.19,
                   f3 = 0.05,
                   f4 = 0,
                   hiv.test.rate = (0.01325/7) * time.unit,
                   test.window.int = 21 / time.unit,
                   art_start = art_start,
                   tx.init.prob = (0.092/7) * time.unit,
                   tx.halt.prob = (0.0102/7) * time.unit,
                   tx.reinit.prob = (0.00066/7) * time.unit,
                   trans.r = 0.05,
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
                   a1.rate = 0.00006 * time.unit,
                   a2.rate = 0.00006 * time.unit,
                   arrival.age = 18,
                   m12.rate = 0,
                   m21.rate = 0)

init <- init.net()

#in weeks
#nsteps = years * 52
# in days
nsteps = years * 365


control <- control.net(type = NULL, nsteps = nsteps, start = 1, nsims = 1,
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


sim <- readRDS("results_sim.RDS")
# Plot epidemiological quantities of interest ----
simdf <- as.data.frame(sim)
plot(simdf$time, simdf$a1.flow)
plot(simdf$time, simdf$a2.flow)
plot(simdf$time, simdf$i.num.pop1/simdf$num.pop1)
plot(simdf$time, simdf$i.num.pop2/simdf$num.pop2)
plot(simdf$time, simdf$i.num.pop1)
plot(simdf$time, simdf$i.num.pop2)

plot(sim, y = c("i.num.pop1", "i.num.pop2"), qnts = 1, legend = TRUE)

plot(sim, y = c("dall_pop1.flow", "dall_pop2.flow"), qnts = 1, legend = TRUE)
plot(sim, y = c("daids_pop1.flow", "daids_pop2.flow"), qnts = 1, legend = TRUE)
plot(sim, y = c("dhiv_pop1.flow", "dhiv_pop2.flow"), qnts = 1, legend = TRUE)

plot(sim, y = c("a1.flow", "a2.flow"), qnts = 1, legend = TRUE)
plot(sim, y = "a1.flow", qnts = 1, legend = TRUE)
plot(sim, y = "a2.flow", qnts = 1, legend = TRUE)
plot(sim, y = c("i.num.pop1", "s.num.pop1"), qnts = 1, legend = TRUE)
plot(sim, y = c("i.num.pop2", "s.num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("i.num.pop1", "s.num.pop1","i.num.pop2", "s.num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = "i.num.pop1", qnts = 1, legend = TRUE)
plot(sim, y = "i.num.pop2", qnts = 1, legend = TRUE)
plot(sim, y = c("num.pop1", "num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("s.num.pop2", "i.num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("s.num.pop1", "i.num.pop1"), qnts = 1, legend = TRUE)
plot(sim, y = c("i.prev.pop1", "i.prev.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = "i.prev.pop2", qnts = 1, legend = TRUE)
plot(sim, y = "i.prev.pop1", qnts = 1, legend = TRUE)
plot(sim, y = c("nArrivals_mig1", "nArrivals_mig2"), qnts = 1, legend = TRUE)
plot(sim, y = "nArrivals_mig1", qnts = 1, legend = TRUE)
plot(sim, y = "nArrivals_mig2", qnts = 1, legend = TRUE)

plot(sim, y = c("incid.pop1", "incid.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("i.num.pop1", "i.num.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = "tot.neg.tests", qnts = 1, legend = TRUE)
plot(sim, y = "newDx", qnts = 1, legend = TRUE)
plot(sim, y= c("hstage0.pop1", "hstage1.pop1", "hstage2.pop1", "hstage3.pop1", "hstage.aids.pop1"), qnts = 1, legend = TRUE)
plot(sim, y= c("hstage0.pop2", "hstage1.pop2", "hstage2.pop2", "hstage3.pop2", "hstage.aids.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("hstage0.pop1", "hstage0.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("hstage1.pop1", "hstage1.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("hstage2.pop1", "hstage2.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("hstage3.pop1", "hstage3.pop2"), qnts = 1, legend = TRUE)
plot(sim, y = c("hstage.aids.pop1", "hstage.aids.pop2"), qnts = 1, legend = TRUE)



# ndtv plots ----
# this plot will only work if setting in control tergmLite = TRUE
# note that creating this plot can be computationally intensive
# so best is to use a small size population or check
# the ndtv package for option to save larger network sizes.
library(ndtv)
slice.par<-list(start=1,end=31,interval=1,
                aggregate.dur=1,rule="earliest")


sim_ani<-compute.animation(sim$network$sim1[[1]],
                           default.dist=3,
                           slice.par=slice.par,
                           animation.mode='MDSJ',
                           verbose=FALSE)


render.d3movie(sim_ani,vertex.col="global_track",
               edge.col="black", edge.lwd = 2,
               displaylabels=TRUE,label.cex=1,
               vertex.cex = 4,
               label.col="blue", verbose=FALSE,
               main='Simulation interactions',
               output.mode = 'htmlWidget')

ani_global <- render.animation(sim_ani,vertex.col="global_track",
                             edge.col="black", edge.lwd = 2,
                             displaylabels=TRUE,label.cex=1,
                             vertex.cex = 4,
                             label.col="blue", verbose=FALSE,
                             main='Simulation interactions',
                             output.mode = 'htmlWidget')


saveGIF(render.animation(sim_ani,vertex.col="global_track",
                         edge.col="black", edge.lwd = 2,
                         displaylabels=TRUE,label.cex=1,
                         vertex.cex = 4,
                         label.col="blue", verbose=FALSE,
                         main='Simulation interactions',
                         output.mode = 'htmlWidget'), "ani_global.gif")

