# plot results of ECDC incidence
library(ggplot2)
library(reshape2)

model1 <- readRDS("../../HIV_platform/Results_modelling1/HIVModelMainFit_20211022_110848.rds")
model2 <- readRDS("../../HIV_platform/Results_modelling2/HIVModelMainFit_20211027_103512.rds")
model3 <- readRDS("../../HIV_platform/Results_modelling3/HIVModelMainFit_20211027_110750.rds")



incidence <- data.frame(year = model1$Year,
                        model1 = model1$N_Inf_M,
                        model2 = model2$N_Inf_M,
                        model3 = model3$N_Inf_M)

incidence_melt <- melt(incidence, id.vars = "year")


ggplot(data = incidence_melt, aes(x = year, y = value, colour = variable)) +
  geom_line() +
  theme_bw()
