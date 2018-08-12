# Title     : p-values —> q-values
# Objective : FDR correction using qvalue package
# Created by: vvi
# Created on: 16.04.18

setwd("~/Science/eQTL_analysis/data/pQTLs/temp/")
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
module_name = args[1]
pvalues_path = sprintf("./pvalues_%s.csv", module_name)
pobj = read.table(pvalues_path, stringsAsFactors=FALSE, col.names=c("pvalues"))
pvalues = as.numeric(pobj$pvalues[- 1])
qobj = try(qvalue(pvalues));
if (class(qobj) == "try-error") {
    qvalues = p.adjust(pvalues, method = "BH")
} else {
    qvalues = qobj$qvalues
}
qvalues_path = sprintf("./qvalues_%s.csv", module_name)
write.table(
qvalues, file = qvalues_path,
col.names = c("qvalue"), row.names = FALSE
)
