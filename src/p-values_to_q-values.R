#!/usr/bin/Rscript

# Title     : p-values —> q-values
# Objective : FDR correction using qvalue package
# Created by: vvi
# Created on: 16.04.18

setwd("~/Science/eQTL_analysis/data/pQTLs/")
library(qvalue)
pobj = read.table("./pvalues.csv", stringsAsFactors=FALSE, col.names=c("pvalues"))
pvalues = as.numeric(pobj$pvalues[- 1])
qobj = try(qvalue(pvalues));
if (class(qobj) == "try-error") {
    qvalues = p.adjust(pvalues, method = "BH")
} else {
    qvalues = qobj$qvalues
}

write.table(
qvalues, file = "./qvalues.csv",
col.names = c("qvalue"), row.names = FALSE
)
