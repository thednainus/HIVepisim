#total number of jobs to add in the DAG file
#this generates a dag file to simulate 30 replicates using the same combination
#of parameter values
n = 30
line_number <- 1067

#lines <- as.list(rep(1:n))
lines <- as.list(rep(line_number, n))
jobs <- as.list(rep(1:n))
dagfilename <- "seed_30jobs.dag"

#part1 <- unlist(lapply(lines, function(x) paste("JOB SIM", x, " hiv_epimodel_simulation2.sub", sep = "")))
#part2 <- unlist(lapply(lines, function(x) paste("VARS SIM", x, " line_number=", "\"", x, "\"",
#                                                " seed=", "\"", floor(runif(1, min=10000, max=100001)), "\"",
#                                                " dirname=\"params_", x, "/rep_1\"", sep = "")))

part1 <- unlist(lapply(jobs, function(x) paste("JOB SIM", x, " hiv_epimodel_simulation2.sub", sep = "")))
part2 <- unlist(lapply(jobs, function(x) paste("VARS SIM", x, " line_number=", "\"", line_number, "\"",
                                               " seed=", "\"", floor(runif(1, min=10000, max=100001)), "\"",
                                               paste(" dirname=\"params_", line_number, sep = ""), "/rep_", x, "\"", sep = "")))

#part1 <- unlist(lapply(jobs, function(x) paste("JOB SIM", x, " 500mig_1067_deepseq.sub", sep = "")))
#part2 <- unlist(lapply(jobs, function(x) paste("VARS SIM", x, " line_number=", "\"", x, "\"",
#                                               " seed=", "\"", floor(runif(1, min=10000, max=100001)), "\"",
#                                               paste(" dirname=\"params_", line_number, sep = ""), "/rep_", x, "\"", sep = "")))


write.table(part1, file = dagfilename, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

write.table(part2, file = dagfilename, quote = FALSE,
            row.names = FALSE, col.names = FALSE,
            append =  TRUE)
