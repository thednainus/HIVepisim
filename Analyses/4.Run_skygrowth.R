# run skygrowth on true and estimated trees
library(skygrowth)
library(ape)
library(ggplot2)

## skygrowth ----

#read tree



tree <- read.tree("treedater_tree.tre")
class(region_only_dated_tree) <- "phylo"

set.seed(12345)
mcmcfit <- skygrowth.mcmc(tree)
mcmcfit$time <- mcmcfit$time + 2021

set.seed(12345)
mcmcfit2 <- skygrowth.mcmc(region_only_dated_tree)
mcmcfit2$time <- mcmcfit2$time + 2021

set.seed(12345)
mcmcfit3 <- skygrowth.mcmc(tree2)
mcmcfit3$time <- mcmcfit3$time + 2021
