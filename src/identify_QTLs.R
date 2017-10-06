# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

identify_QTLs <- function(genotypes_filename, expression_filename, output_filename) {
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
  
  # Only associations significant at this level will be saved
  pvOutputThreshold = 1e-3;
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = '\t';      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(full_genotypes_filename);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = '\t';      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(full_expression_filename);
  
  ## Run the analysis
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    output_file_name = full_output_filename,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = TRUE);
  
  ## Results:
  
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  cat('Detected QTLs:', '\n');
  
  ## Plot the histogram of all p-values
  
  plot(me)
}

identify_QTLs('rna_genotypes.csv', 'rna_expression_avg.csv', "eQTLs.csv")
identify_QTLs('protein_genotypes.csv', 'protein_expression_avg.csv', "pQTLs.csv")