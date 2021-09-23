lines <- as.list(rep(1:1000))
dagfilename <- "largepop_sim.dag"

part1 <- unlist(lapply(lines, function(x) paste("JOB SIM", x, " hiv_epimodel_simulation.sub", sep = "")))
part2 <- unlist(lapply(lines, function(x) paste("VARS SIM1 line_numer=", x, sep = "")))


write.table(c(part1, part2), file = dagfilename, quote = FALSE,
            row.names = FALSE, col.names = FALSE)
