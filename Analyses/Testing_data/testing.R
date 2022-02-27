#hiv.test.rate:data was previously fit using a fixed value of 0.045
total_testers <- data.frame(year=c(1985:2001),
                            total_people_tested=c(2118,2888, 13012, 11451, 9682,
                                                  17428, 25565, 28826, 19111, 15456,
                                                  17563, 17208, 16730, 15836, 17351,
                                                  17156, 17077))
hiv.test.rate.df <- readRDS(system.file("data/probability_msm_tested_per_day.RDS",
                                        package = "HIVepisim"))

total_testers1990_to_2001 <- total_testers[6:17,]
total_testers1990_to_2001["total_males_tested"] <- hiv.test.rate.df$total_males_tested
total_testers1990_to_2001["total_msm_tested"] <- hiv.test.rate.df$total_msm_tested
total_testers1990_to_2001["perc_males_tested"] <- total_testers1990_to_2001$total_males_tested/total_testers1990_to_2001$total_people_tested
total_testers1990_to_2001["perc_msm_tested"] <- total_testers1990_to_2001$total_msm_tested/total_testers1990_to_2001$total_people_tested
#on average based on data fro 1990 to 2001 I have 43.6% are males getting tested
average_males_tested <- mean(total_testers1990_to_2001$perc_males_tested)
#on average based on data fro 1990 to 2001 I have 11.7% are males getting tested
average_msm_tested <- mean(total_testers1990_to_2001$perc_msm_tested)

total_testers1985_1989 <- total_testers[1:5,]
total_testers1985_1989["total_males_tested"] <- average_males_tested * total_testers1985_1989$total_people_tested
total_testers1985_1989["total_msm_tested"] <- average_msm_tested * total_testers1985_1989$total_people_tested
total_testers1985_1989["perc_msm_tested"] <- total_testers1985_1989$total_msm_tested/total_testers1985_1989$total_males_tested


#based on this answer:
#https://math.stackexchange.com/questions/357242/calculating-probabilities-over-different-time-intervals
total_testers1985_1989["perc_per_day"] <- 1 - ((1-total_testers1985_1989$perc_msm_tested)^(1/365))


# this is the data Sanjay Mehta shared with us on 26 January 2022

library(reshape)

testing <- read.csv("~/Desktop/Imperial/newHIVproject-01Aug2020/SanDiego_testing_data/SDPHTESTING.UPDATED.HHSA2_July19Update.csv")

total_tests <- apply(testing[,3:21], 2, sum, na.rm = TRUE)
names(total_tests) <- c(1998:2016)

total_tests <- data.frame(year = names(total_tests), total = unname(total_tests))

