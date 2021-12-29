lines <- as.list(rep(1:1000))
dagfilename <- "largepop_sim.dag"

part1 <- unlist(lapply(lines, function(x) paste("JOB SIM", x, " hiv_epimodel_simulation.sub", sep = "")))
part2 <- unlist(lapply(lines, function(x) paste("VARS SIM", x, " line_number=", "\"", x, "\"",
                                                " dirname=\"params_", x, "/rep_1\"", sep = "")))


write.table(part1, file = dagfilename, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

write.table(part2, file = dagfilename, quote = FALSE,
            row.names = FALSE, col.names = FALSE,
            append =  TRUE)
