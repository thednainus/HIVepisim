#hiv.test.rate:data was previously fit using a fixed value of 0.045
#create table for probability that an MSN will be tested
#data from HIV Trends and the Status of High Risk Groups in San Diego County, 1985-2001
#in https://www.sandiegocounty.gov/content/sdc/hhsa/programs/phs/hiv_aids_epidemiology_unit/archives.html

test_msm_data <- data.frame(year = c(1990:2001),
                            total_msm_tested = c(2831, 3901, 3163, 2510, 1932,
                                                 1908, 1354, 1384, 1664, 1752,
                                                 1807, 2488),
                            total_bisexual_tested =c(589, 882, 670, 514, 549, 701,
                                                     781, 777, 456, 532, 599, 482),
                            total_testers = c(11078, 18203, 19189, 14501, 11946,
                                              12589, 12710, 12107, 10274, 12262,
                                              11856, 10772))
test_msm_data["total_msm"] <- test_msm_data$total_msm_tested + test_msm_data$total_bisexual_tested

test_msm_data["perc_msm_tested"] <- test_msm_data$total_msm/test_msm_data$total_testers


total_testers <- data.frame(year=c(1985:2001),
                            total_people_tested=c(2118,2888, 13012, 11451, 9682,
                                                  17428, 25565, 28826, 19111, 15456,
                                                  17563, 17208, 16730, 15836, 17351,
                                                  17156, 17077),
                            total_anonymous_testers = c(rep(NA, 5),11078, 18203, 19189, 14501,
                                                        11946, 12589, 12710, 12107,
                                                        10274, 12262, 11856, 10772))

#average of proportions
average_in_anonymous <- sum(total_testers$total_anonymous_testers[6:17])/sum(total_testers$total_people_tested[6:17])
total_testers$total_anonymous_testers[1:6] <- total_testers$total_people_tested[1:6] *average_in_anonymous

#average number of MSM (percentage of MSm tested from)
average_msm <- sum(test_msm_data$total_msm)/sum(test_msm_data$total_testers)

#take this average for 1985 to 1989
new_data <- total_testers[1:6,]
new_data$total_msm <- new_data$total_anonymous_testers * average_msm
new_data["perc_msm_tested"] <- new_data$total_msm/new_data$total_anonymous_testers

test_msm_data1989 <- data.frame(year=1989, total_msm_tested = NA,
                                total_bisexual_tested = NA, total_testers = 6767.607,
                                total_msm  = 1470.776, perc_msm_tested = 0.2173259)
test_msm_data <- rbind(test_msm_data1989, test_msm_data)

#based on this answer:
#https://math.stackexchange.com/questions/357242/calculating-probabilities-over-different-time-intervals
test_msm_data["perc_per_day"] <- 1 - ((1-test_msm_data$perc_msm_tested)^(1/365))
saveRDS(test_msm_data, "inst/data/probability_msm_tested_per_day2.RDS")





# this is the data Sanjay Mehta shared with us on 26 January 2022

library(reshape)

testing <- read.csv("~/Desktop/Imperial/newHIVproject-01Aug2020/SanDiego_testing_data/SDPHTESTING.UPDATED.HHSA2_July19Update.csv")

total_tests <- apply(testing[,3:21], 2, sum, na.rm = TRUE)
names(total_tests) <- c(1998:2016)

total_tests <- data.frame(year = names(total_tests), total = unname(total_tests))

