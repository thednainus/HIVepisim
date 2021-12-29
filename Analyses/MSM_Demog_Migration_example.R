library(EpiModel)
#library(HIVepisim)
#library(EpiModelHPC)

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

dr_vec <- rep(dr_pp_pw, times = age_spans)
data.frame(ages, dr_vec)

# Plot death rates
#par(mar = c(3,3,2,1), mgp = c(2,1,0), mfrow = c(1,1))
#plot(ages, dr_vec, type = "o", xlab = "age", ylab = "Mortality Rate")

# Initialize network
#n_pop1 = 50000
#n_pop2 = 300000
n_pop1 = 1000
n_pop2 = 2000

# total number of individuals in network
n_total = n_pop1 + n_pop2

#initialize network
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
n_inf_pop1 <- 4
init.Infected1 <- sample(1:n_pop1, n_inf_pop1)
diagStatusVec1[init.Infected1] <- 1


diagStatusVec2 <- rep(0, n_pop2)
n_inf_pop2 <- 10
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

# When forcing degree = 1 only (no concurrency)
#formation <- ~edges + nodemix("origin", levels2 = c(1, 2)) + degree(0, by = "origin") +
#  concurrent(by="origin")

#simplest formation model
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

# Dissolution model
durs <- c(23.7, 23.7, 2)
coef.diss <- dissolution_coefs(~offset(edges) +
                               offset(nodemix("origin", levels2 = c(1, 2))),
                               duration = durs, d.rate = mean(dr_vec))
coef.diss

# Fit the model
# # Fit the TERGM
# to simulate large networks
# est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE,
#              set.control.ergm = control.ergm(MCMC.burnin=1e7, MCMC.interval=1e7,
#                                              MCMC.samplesize=20000,
#                                              init.MPLE.samplesize = 1e8,
#                                              SAN.control = control.san(SAN.nsteps = 1e8)
#                                              ))

est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE)


# Model diagnostics
# Simulate time series to examine timed edgelist
# to simulate large networks
#dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE,
#            set.control.ergm = control.simulate.ergm(MCMC.init.maxchanges = 1e8,
#                                                     MCMC.burnin.min= 1.5e5,
#                                                     MCMC.burnin.max =1.5e5))


dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE)

#dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE,
#            nwstats.formula = ~edges + nodemix("origin", levels2 = c(1, 2)) + degree(0:3, by = "origin"))
#dx
#plot(dx)
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




