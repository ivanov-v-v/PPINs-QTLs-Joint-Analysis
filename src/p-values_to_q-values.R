#!/usr/bin/Rscript

# Title     : p-values —> q-values
# Objective : FDR correction using qvalue package
# Created by: vvi
# Created on: 16.04.18

setwd("~/Science/eQTL_analysis/data/pQTLs/")
library(qvalue)
pobj = read.table("./pvalues.csv", stringsAsFactors=FALSE, col.names=c("pvalues"))
qobj = qvalue(as.numeric(pobj$pvalues[-1]))
write.table(
    qobj$qvalues, file="./qvalues.csv",
    col.names=c("qvalue"), row.names=FALSE
)
