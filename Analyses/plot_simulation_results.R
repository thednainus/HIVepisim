sim1 <- readRDS("~/Desktop/tmp2/osg_new_scrips/results/params_2575/rep_1/results/results_sim.RDS")
nodes1 <- read.csv("~/Desktop/tmp2/osg_new_scrips/results/params_2575/rep_1/results/total_nodes.csv")
nodes1 <- rbind(c(1,60000,20000), nodes1)
sim2 <- readRDS("~/Desktop/tmp2/osg_new_scrips/results/params_6284/rep_1/results/results_sim.RDS")
nodes2 <- read.csv("~/Desktop/tmp2/osg_new_scrips/results/params_6284/rep_1/results/total_nodes.csv")
nodes2 <- rbind(c(1,60000,20000), nodes2)


sim1$stats$nwstats$sim1[[1]] <- cbind(sim1$stats$nwstats$sim1[[1]] , global_meandeg = sim1$stats$nwstats$sim1[[1]][,5]/nodes1$global)
sim1$stats$nwstats$sim1[[1]] <- cbind(sim1$stats$nwstats$sim1[[1]] , region_meandeg = sim1$stats$nwstats$sim1[[1]][,6]/nodes1$region)

sim1$stats$nwstats$sim1[[1]] <- cbind(sim1$stats$nwstats$sim1[[1]] , global_meandeg2 = (sim1$stats$nwstats$sim1[[1]][,2]*2)/nodes1$global)
sim1$stats$nwstats$sim1[[1]] <- cbind(sim1$stats$nwstats$sim1[[1]] , region_meandeg2 = (sim1$stats$nwstats$sim1[[1]][,4]*2)/nodes1$region)


sim2$stats$nwstats$sim1[[1]] <- cbind(sim2$stats$nwstats$sim1[[1]] , global_meandeg = sim2$stats$nwstats$sim1[[1]][,5]/nodes2$global)
sim2$stats$nwstats$sim1[[1]] <- cbind(sim2$stats$nwstats$sim1[[1]] , region_meandeg = sim2$stats$nwstats$sim1[[1]][,6]/nodes2$region)

plot(sim1, type = "formation", plots.joined = FALSE, stats = c("global_meandeg", "region_meandeg"))
plot(sim1, type = "formation", plots.joined = TRUE, stats = c("global_meandeg", "global_meandeg2"))
plot(sim1, type = "formation", plots.joined = TRUE, stats ="meandeg")
plot(sim1, type = "formation", plots.joined = TRUE, stats = c("region_meandeg", "region_meandeg2"))
plot(sim2, type = "formation", plots.joined = FALSE, stats = c("global_meandeg", "region_meandeg"))



plot(sim1, y=c("dall_pop1.flow", "dall_pop2.flow",
              "nArrivals_mig1", "nArrivals_mig2",
              "a1.flow", "a2.flow"), legend = TRUE)

plot(sim1, y=c("i.num.pop1", "s.num.pop1"), legend = TRUE)
plot(sim1, y=c("i.num.pop2", "s.num.pop2"), legend = TRUE)
plot(sim, y=c("i.num.pop1", "s.num.pop1", "i.num.pop2", "s.num.pop2"), legend = TRUE)
