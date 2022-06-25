terminalMTFixed <- read.csv("../../Result/Simulation/MTFixed/out_terminal.csv")
terminalMTVariable <- read.csv("../../Result/Simulation/MTVariable/out_terminal.csv")

outFile <- "../../Result/Analysis/StatisticalTest/wilcoxon_ranksum_test.txt"

sink( outFile)
cat( "Wilcoxon Rank Sum Test : Movement of MT Fixed and MT Variable")
print(wilcox.test( terminalMTFixed$Nuc_distance.per. , terminalMTVariable$Nuc_distance.per. ))
cat("Wilcoxon Rank Sum Test : Angle of MT Fixed and MT Variable")
print(wilcox.test( terminalMTFixed$angel.degree. , terminalMTVariable$angel.degree. ))
sink()