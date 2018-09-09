# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Be sure to use an up to date version of R and Matrix eQTL.

identify_QTLs <-
  function(genotypes_filename,
           expression_filename,
           output_filename,
           covariates_filename = NULL,
           P_THRESHOLD = 0.05) {
    library(MatrixEQTL)
    
    ## Location of the package with the data files.
    base.dir = "~/Science/eQTL_analysis/data/"
    
    
    ## Settings
    
    # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
    useModel = modelANOVA
    # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
    
    # Genotype file name
    full_genotypes_filename = paste(base.dir, genotypes_filename, sep = "")
    
    # Gene expression file name
    full_expression_filename = paste(base.dir, expression_filename, sep = "")
    
    # Covariates file name
    full_covariates_filename = paste(base.dir, covariates_filename, sep="")
    
    # Output file name
    full_output_filename = paste(base.dir, output_filename, sep = "")
    
    if (!file.exists(full_output_filename)) {
      file.create(full_output_filename)
    }
    
    # Only associations significant at this level will be saved
    pvOutputThreshold = P_THRESHOLD
    
    
    ## Load genotype data
    
    snps = SlicedData$new()
    
    snps$fileDelimiter = '\t'
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 5000
    snps$LoadFile(full_genotypes_filename)
    
    
    ## Load gene expression data
    
    gene = SlicedData$new()
    
    gene$fileDelimiter = '\t'
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 5000
    gene$LoadFile(full_expression_filename)
    
    
    ## Load covariates data
    
    covariates = SlicedData$new()
    # default argument (from the man page)
    if (!is.null(covariates_filename)) {
      covariates$fileDelimiter = '\t'
      covariates$fileOmitCharacters = "NA"
      covariates$fileSkipRows = 1
      covariates$fileSkipColumns = 1
      covariates$fileSliceSize = 5000
      covariates$LoadFile(full_covariates_filename)
      
    }
    
    ## Run the analysis
    
    me = Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      cvrt = covariates,
      output_file_name = full_output_filename,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      verbose = TRUE,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )
    
    
    ## Results:
    
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
  }

# args <- commandArgs(TRUE)
library(qvalue)
setwd("~/Science/eQTL_analysis/data/")

genotypes_filename = 'eQTLs/2018/genotypes_filtered_MatrixEQTL.csv'
expression_filename = 'eQTLs/2018/expression.csv'
covariates_filename = 'eQTLs/2018/covariates_MatrixEQTL.csv'
output_filename = 'eQTLs/2018/qtls_MatrixEQTL_with_covariates_ANOVA.csv'

identify_QTLs(
  genotypes_filename = genotypes_filename,
  expression_filename = expression_filename,
  covariates_filename = covariates_filename,
  output_filename = output_filename,
  P_THRESHOLD = 0.05
)

QTLs_df = na.omit(
  read.table(
    output_filename,
    sep = '\t',
    stringsAsFactors = FALSE,
    header = TRUE,
    check.names = FALSE
  )
)
QTLs_df$q.value = qvalue(
  p=as.numeric(QTLs_df$p.value),
  lambda=seq(from=0, to=max(QTLs_df$p.value), 0.05)
)$qvalues
names(QTLs_df)[names(QTLs_df) == "q.value"] = "q_value"
names(QTLs_df)[names(QTLs_df) == "p.value"] = "p_value"
QTLs_df = QTLs_df[, c("SNP", "gene", "p_value", "q_value")]
write.table(
  QTLs_df[QTLs_df$q_value <= 0.05, ],
  output_filename,
  sep = '\t',
  row.names = FALSE
)
gc()
