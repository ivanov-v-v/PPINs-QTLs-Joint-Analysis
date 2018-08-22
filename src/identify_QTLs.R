# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Be sure to use an up to date version of R and Matrix eQTL.

identify_QTLs <- function(genotypes_filename, expression_filename, output_filename, P_THRESHOLD = 1) {
  library(MatrixEQTL)
  
  ## Location of the package with the data files.
  base.dir = "~/Science/eQTL_analysis";
  
  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR;# modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  full_genotypes_filename = paste(base.dir, paste("/data/", genotypes_filename, sep=""), sep="");
  
  # Gene expression file name
  full_expression_filename = paste(base.dir, paste("/data/", expression_filename, sep=""), sep="");
  
  # Output file name
  full_output_filename = paste(base.dir, paste("/data/", output_filename, sep=""), sep="");
  if (!file.exists(full_output_filename)) {
    file.create(full_output_filename)
  }
  
  # Only associations significant at this level will be saved
  pvOutputThreshold = P_THRESHOLD;
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = '\t';      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 5000;      # read file in slices of 2,000 rows
  snps$LoadFile(full_genotypes_filename);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = '\t';      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 5000;      # read file in slices of 5,000 rows
  gene$LoadFile(full_expression_filename);
  
  ## Run the analysis
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    output_file_name = full_output_filename,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    verbose = TRUE,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  );
  
  ## Results:
  
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  
  ## Plot the histogram of all p-values
}

# args <- commandArgs(TRUE)

identify_QTLs(
'eQTLs/genotypes_2018_filtered.csv',
'eQTLs/expression_2018.csv',
'eQTLs/qtls_interpolated_MatrixEQTL.csv',
 0.05
)

library(qvalue)

QTLs_df = na.omit(read.table("~/Science/eQTL_analysis/data/eQTLs/qtls_interpolated_MatrixEQTL.csv", sep = '\t', stringsAsFactors = FALSE, header = TRUE))
QTLs_df$q.value = qvalue(
  p=as.numeric(QTLs_df$p.value),
  lambda=seq(from=0, to=max(QTLs_df$p.value), 0.05)
)$qvalues
#QTLs_df = QTLs_df[, c(1, 2, 5 : 7)]
write.table(QTLs_df[QTLs_df$q.value <= 0.05,], "~/Science/eQTL_analysis/data/eQTLs/qtls_interpolated_MatrixEQTL.csv", sep = '\t', row.names = FALSE)
