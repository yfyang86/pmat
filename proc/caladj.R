#!/usr/local/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1], header = F)
colnames(data) <- c("chr", "start", "stop", "abs.methyl.diff", "CpGs", "pFN", "5mC.A", "5mC.B")
data$FDR <- p.adjust(data$pFN, method = "BH")
write.table(data, file = "test.mr.DMR.fdr", row.names= F, sep = "\t", quote = F)

