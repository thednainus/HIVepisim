#sample parameters to simulate MSM epidemics
library(lhs)
library(HIVepisim)
sampler6 <- function(n = 100, paramdf){

  ## number of dimension (the number of parameters that we will estimate)
  d = length(names(paramdf))
  params <- randomLHS(n, d)
  params[,1] <- qunif(params[,1], min = paramdf[1,1], max = paramdf[2,1])
  params[,2] <- qunif(params[,2], min = paramdf[1,2], max = paramdf[2,2])
  params[,3] <- qunif(params[,3], min = paramdf[1,3], max = paramdf[2,3])
  params[,4] <- qunif(params[,4], min = paramdf[1,4], max = paramdf[2,4])
  params[,5] <- qunif(params[,5], min = paramdf[1,5], max = paramdf[2,5])
  params[,6] <- qunif(params[,6], min = paramdf[1,6], max = paramdf[2,6])



  params <- data.frame(params)
  colnames(params) <- names(paramdf)


  return (params)

}


# parameters to sample:
# We use 3 different sclar number to try to fit the incidence curve
# that has a sharp decrease around 1995
# 1. Scalar to increase act.rate so as to obtain a
# number of infected individuals equivalent to real data

# 2. Scalar to increase act.rate so as to obtain a
# number of infected individuals equivalent to real data

# 3. Scalar to increase act.rate so as to obtain a
# number of infected individuals equivalent to real data

# 4. Probability of initiating treatment per day: # 63% per year is 0.0027 per day
# 99% per year is 0.01253765 per day.

# 5. Initial size of infected population for the MSM in San Diego

# for OSG
params6dim <- sampler6(n = 10000,
                       paramdf = data.frame(scalar1 =  c(3000, 6000),
                                            scalar2 =  c(3000, 7000),
                                            scalar3 =  c(1500, 4000),
                                            tx.init.prob.param = c(0.0001, 0.01253765),
                                            init_pop1_param = c(10, 300),
                                            scalar4 =  c(1, 3)))


#for imperial college cluster using transmission rate per year (intead of a
# contant rate)
params6dim <- sampler6(n = 10000,
                       paramdf = data.frame(scalar1 =  c(800, 2000),
                                            scalar2 =  c(20000, 30000),
                                            scalar3 =  c(10000, 25000),
                                            tx.init.prob.param = c(0.0001, 0.01253765),
                                            init_pop1_param = c(10, 300),
                                            scalar4 =  c(1, 3)))

# for imperial college cluster but using a certain paramter values that seem
# to be closer to the diagnosis incidence (some preminaly analysis in OSG)
params6dim <- sampler6(n = 5000,
                       paramdf = data.frame(scalar1 =  c(3000, 3500),
                                            scalar2 =  c(4000, 5300),
                                            scalar3 =  c(3500, 5000),
                                            tx.init.prob.param = c(0.0001, 0.01253765),
                                            init_pop1_param = c(40, 300),
                                            scalar4 =  c(1, 2)))

saveRDS(params6dim, "inst/data/params6dim_21Feb2022_xsede_5000jobs.RDS")

