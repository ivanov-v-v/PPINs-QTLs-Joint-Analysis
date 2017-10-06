All the data was originally used in the series of papers from Kruglyak's lab in Princeton:
- [Genetic Dissection of Transcriptional Regulation in Budding Yeast](http://science.sciencemag.org/content/296/5568/752)
- [The landscape of genetic complexity across 5,700 gene expression traits in yeast, 2005](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC547855/)

## rna_expression_full.csv
Processed and cleaned up data used in the original paper.
- Column 0: indices.
- Column 1: gene names.
- Columns 2-7: BY parent.
- Columns 8-19: RM parent.
- Columns 20-132: progeny expression values.

Compiled from GDS1115 and GDS1116 data in the following way:
* For each gene, for each strain, average the measurements obtained from two microarrays
	- If both values are available, average them
	- If only one value is available, take that value
	- If both values are not present, write NaN
* Drop the rows without annotation (either "empty" or "blank")
* Replace GSM labels with exact strain names. 
* Drop the columns corresponding to non-genotyped strains
* Sort all the progeny strains accordingly to their name, tokenized into 3-tuples of the form (int, int, str)

## rna_expression_avg.csv
The same as rna_expression_full.csv, with parent strains averaged.
- Column 0: indices
- Column 1: gene names.
- Column 2-3: averaged expression of parent strains, BY and RM
- Column 4-128: progeny expression values

## protein_expression_avg.csv

## genotypes_full.csv:
Genome annotation with 2820 genetic markers for budding yeast strains BY and RM and their progeny, 124 genotypes in total. Covers all the strains mentioned in the original Brem paper and also some extra ones.

## rna_genotypes

## protein_genotypes.csv

## eQTLs.csv

## pQTLs.csv
