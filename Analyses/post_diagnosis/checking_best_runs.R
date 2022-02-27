# simulatioon
library(EpiModel)
library(HIVepisim)
library(HIVepisimAnalysis)
library(lubridate)
library(ggplot2)
library(reshape2)

sim1 <- readRDS("/Users/user/Desktop/tmp2/osg_new_scrips/results_narrow_parameters/params_2348/rep_1/results/results_sim.RDS")
sim1_df <- as.data.frame(sim1)


sim2 <- readRDS("/Users/user/Desktop/tmp2/osg_new_scrips/results_narrow_parameters/params_1067/rep_1/results/results_sim.RDS")
sim2_df <- as.data.frame(sim2)

#arrivals in pop1
arrival_sim1_pop1 <- data.frame(time = sim1_df$time, a1.flow = sim1_df$a1.flow)
arrival_sim2_pop1 <- data.frame(time = sim2_df$time, a1.flow = sim2_df$a1.flow)


#arrivals in pop2
arrival_sim1_pop2 <- data.frame(time = sim1_df$time, a2.flow = sim1_df$a2.flow)
arrival_sim2_pop2 <- data.frame(time = sim2_df$time, a2.flow = sim2_df$a2.flow)

#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

#convert days to years
arrival_sim1_pop1["time_years"] <- days2years(arrival_sim1_pop1$time,
                                              init_date = init_sim_date)
arrival_sim1_pop1["year"] <- unlist(lapply(arrival_sim1_pop1$time_years,
                                           function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
arrival_sim1_pop1$year <- as.factor(arrival_sim1_pop1$year)


arrival_sim2_pop1["time_years"] <- days2years(arrival_sim2_pop1$time, init_date = init_sim_date)
arrival_sim2_pop1["year"] <- unlist(lapply(arrival_sim2_pop1$time_years,
                                           function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
arrival_sim2_pop1$year <- as.factor(arrival_sim2_pop1$year)

arrival_sim1_pop2["time_years"] <- days2years(arrival_sim1_pop2$time, init_date = init_sim_date)
arrival_sim1_pop2["year"] <- unlist(lapply(arrival_sim1_pop2$time_years,
                                           function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
arrival_sim1_pop2$year <- as.factor(arrival_sim1_pop2$year)


arrival_sim2_pop2["time_years"] <- days2years(arrival_sim2_pop2$time, init_date = init_sim_date)
arrival_sim2_pop2["year"] <- unlist(lapply(arrival_sim2_pop2$time_years,
                                           function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
arrival_sim2_pop2$year <- as.factor(arrival_sim2_pop2$year)


#aggregate arrival by year

arrival_sim1_pop1 <- aggregate(a1.flow  ~ year, data = arrival_sim1_pop1[,c(4,2)], FUN=sum)
arrival_sim2_pop1 <- aggregate(a1.flow  ~ year, data = arrival_sim2_pop1[,c(4,2)], FUN=sum)

arrival_pop1 <- data.frame(year = arrival_sim1_pop1$year,
                           sim1_arrival = arrival_sim1_pop1$a1.flow,
                           sim2_arrival = arrival_sim2_pop1$a1.flow)
arrival_pop1$year <- as.character(arrival_pop1$year)
arrival_pop1$year <- as.numeric(arrival_pop1$year)


arrival_sim1_pop2 <- aggregate(a2.flow  ~ year, data = arrival_sim1_pop2[,c(4,2)], FUN=sum)
arrival_sim2_pop2 <- aggregate(a2.flow  ~ year, data = arrival_sim2_pop2[,c(4,2)], FUN=sum)


arrival_pop2 <- data.frame(year = arrival_sim1_pop2$year,
                           sim1_arrival = arrival_sim1_pop2$a2.flow,
                           sim2_arrival = arrival_sim2_pop2$a2.flow)
arrival_pop2$year <- as.character(arrival_pop2$year)
arrival_pop2$year <- as.numeric(arrival_pop2$year)


#melt data
arrival_pop1m <- melt(arrival_pop1, id.vars = c("year"))

quartz()
ggplot(arrival_pop1m, aes(x = year, y = value)) +
  geom_line(aes(colour = variable), size = 1.5) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20)) +
  labs(y = "arrivals pop1")


arrival_pop2m <- melt(arrival_pop2, id.vars = c("year"))

quartz()
ggplot(arrival_pop2m, aes(x = year, y = value)) +
  geom_line(aes(colour = variable), size = 1.5) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20)) +
  labs(y = "arrivals pop2")


plot(sim1, y=c("i.num.pop1", "s.num.pop1"), legend = TRUE)
plot(sim1, y=c("i.num.pop2", "s.num.pop2"), legend = TRUE)
plot(sim1, y=c("i.num.pop1", "s.num.pop1", "i.num.pop2", "s.num.pop2"), legend = TRUE)

plot(sim1, y=c("hstage0.pop1", "hstage1.pop1", "hstage2.pop1",
              "hstage3.pop1", "hstage.aids.pop1"), legend = TRUE)
plot(sim1, y=c("hstage0.pop1", "hstage1.pop1"), legend = TRUE)
plot(sim1, y=c("hstage.aids.pop1"), legend = TRUE)

plot(sim1, y=c("hstage0.pop2", "hstage1.pop2",
              "hstage2.pop2", "hstage3.pop2",
              "hstage.aids.pop2"), legend = TRUE)
