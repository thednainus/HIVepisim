# plot results of ECDC incidence
library(ggplot2)
library(reshape2)

model1 <- readRDS("../../HIV_platform/Results_modelling1/HIVModelMainFit_20211022_110848.rds")
model2 <- readRDS("../../HIV_platform/Results_modelling2/HIVModelMainFit_20211027_103512.rds")
model3 <- readRDS("../../HIV_platform/Results_modelling3/HIVModelMainFit_20211027_110750.rds")
#6knots
#spline every 2 years from 1992 for the diagnosis rate
model4 <- readRDS("../../HIV_platform/Results_modelling4/HIVModelMainFit_20220131_142936.rds")
#4knots
#spline every 2 years from 1992 for the diagnosis rate
model5 <- readRDS("../../HIV_platform/Results_modelling5/HIVModelMainFit_20220131_145158.rds")



incidence <- data.frame(year = model1$Year,
                        model1 = model1$N_Inf_M,
                        model2 = model2$N_Inf_M,
                        model3 = model3$N_Inf_M,
                        model4 = model4$N_Inf_M,
                        model5 = model5$N_Inf_M)


incidence_melt <- melt(incidence, id.vars = "year")


ggplot(data = incidence_melt, aes(x = year, y = value, colour = variable)) +
  geom_line() +
  theme_bw()


incidence <- data.frame(year = model1$Year,
                        model1 = model1$N_Inf_M)

ggplot(data = incidence, aes(x = year, y = model1)) +
  geom_line(size = 1) +
  theme_bw() +
  xlab("Year") +
  ylab("Incidence") +
  theme(text=element_text(size=30))

#ggsave("incidence_plot_ecdc.pdf", useDingbats=FALSE, width=19.3, height=11.1)
ggsave("incidence_plot_ecdc.pdf", useDingbats=FALSE, width=12, height=9)
