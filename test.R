library(EpiModel)
rm(list = ls())

ages <- 0:85

# Rates per 100,000 for age groups: <1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29,
#                                   30-34, 35-39, 40-44, 45-49, 50-54, 55-59,
#                                   60-64, 65-69, 70-74, 75-79, 80-84, 85+
departure_rate <- c(588.45, 24.8, 11.7, 14.55, 47.85, 88.2, 105.65, 127.2,
                    154.3, 206.5, 309.3, 495.1, 736.85, 1051.15, 1483.45,
                    2294.15, 3642.95, 6139.4, 13938.3)


dr_pp_pd <- departure_rate / 1e5 / 365


age_spans <- c(1, 4, rep(5, 16), 1)
dr_vec <- rep(dr_pp_pd, times = age_spans)
data.frame(ages, dr_vec)


n <- 1000
nw <- network_initialize(n)

ageVec <- sample(ages, n, replace = TRUE)
nw <- set_vertex_attribute(nw, "age", ageVec)


statusVec <- rep("s", n)
init.latent <- sample(1:n, 50)
statusVec[init.latent] <- "e"
statusVec

table(statusVec)

nw <- set_vertex_attribute(nw, "status", statusVec)
nw



aging <- function(dat, at) {

  age <- get_attr(dat, "age")
  age <- age + 1/365
  dat <- set_attr(dat, "age", age)

  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))

  return(dat)
}



dfunc <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters
  dep.rates <- get_param(dat, "departure.rates")
  dep.dis.mult <- get_param(dat, "departure.disease.mult")

  ## Query alive
  idsElig <- which(active == 1)
  nElig <- length(idsElig)

  ## Initialize trackers
  nDepts <- 0
  idsDepts <- NULL

  if (nElig > 0) {

    ## Calculate age-specific departure rates for each eligible node ##
    ## Everyone older than 85 gets the final mortality rate
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    drates_of_elig <- dep.rates[whole_ages_of_elig]

    ## Multiply departure rates for diseased persons
    idsElig.inf <- which(status[idsElig] == "i")
    drates_of_elig[idsElig.inf] <- drates_of_elig[idsElig.inf] *
      dep.dis.mult

    ## Simulate departure process
    vecDepts <- which(rbinom(nElig, 1, drates_of_elig) == 1)
    idsDepts <- idsElig[vecDepts]
    nDepts <- length(idsDepts)

    ## Update nodal attributes
    if (nDepts > 0) {
      active[idsDepts] <- 0
      exitTime[idsDepts] <- at
      dat <- set_attr(dat, "active", active)
      dat <- set_attr(dat, "exitTime", exitTime)
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "total.deaths", at, nDepts)

  # covid deaths
  covid.deaths <- length(intersect(idsDepts, which(status == "i")))
  dat <- set_epi(dat, "covid.deaths", at, covid.deaths)

  return(dat)
}



afunc <- function(dat, at) {

  ## Parameters ##
  n <- get_epi(dat, "num", at - 1)
  a.rate <- get_param(dat, "arrival.rate")

  ## Process ##
  nArrivalsExp <- n * a.rate
  nArrivals <- rpois(1, nArrivalsExp)

  # Update attributes
  if (nArrivals > 0) {
    dat <- append_core_attr(dat, at = at, n.new = nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)

    dat <- append_attr(dat, "age", 0, nArrivals)
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}


formation <- ~edges + degree(0) + absdiff("age")


mean_degree <- 2
edges <- mean_degree * (n/2)
avg.abs.age.diff <- 2
isolates <- n * 0.08
absdiff <- edges * avg.abs.age.diff

target.stats <- c(edges, isolates, absdiff)
target.stats

coef.diss <- dissolution_coefs(~offset(edges), 20, mean(dr_vec))
coef.diss


est <- netest(nw, formation, target.stats, coef.diss)


dx <- netdx(est, nsims = 10, ncores = 5, nsteps = 500,
            nwstats.formula = ~edges + absdiff("age") + degree(0:6))

print(dx)
plot(dx)

param <- param.net(inf.prob = 0.0004,
                   act.rate = 1,
                   departure.rates = dr_vec,
                   departure.disease.mult = 100,
                   arrival.rate = 3.223207e-04,
                   ei.rate = 0.05, ir.rate = 0)


init <- init.net()

control <- control.net(type = NULL,
                       nsims = 1,
                       ncores = 1,
                       nsteps = 2000,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       aging.FUN = aging,
                       departures.FUN = dfunc,
                       arrivals.FUN = afunc,
                       resimulate.network = TRUE,
                       tergmLite = TRUE)


sim <- netsim(est, param, init, control)
tm1 <- get_transmat(sim)


el <- cbind(tm1$inf, tm1$sus)
IDPOP <- unique(as.vector(el))

sus_unique <- length(unique(tm1$sus))
sus_all <- length(tm1$sus)

inf_unique <- length(unique(tm1$inf))
inf_all <- length(tm1$inf)
