#!/usr/local/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1], header = T)
## users could change the p-value adjustment method
if (length(args) == 1){
    method = "fdr"
}else{
    if (args[2] %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){
        method = args[2]
    }else{
        method = "fdr"
    }
}
data$FDR <- p.adjust(data$pFN, method = method)
names(data)[ncol(data)] = toupper(method)
write.table(data, file = paste("test.mr.DMR.", method, sep = '') , 
            row.names= F, sep = "\t", quote = F)
