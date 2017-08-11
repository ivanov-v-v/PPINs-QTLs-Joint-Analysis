## brem2005_strains.csv
Names of strains used in a series of papers from Kruglyak's lab:
- [Genetic Dissection of Transcriptional Regulation in Budding Yeast](http://science.sciencemag.org/content/296/5568/752)
- [The landscape of genetic complexity across 5,700 gene expression traits in yeast, 2005](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC547855/)

## jbloom_strains_genotyped.tab:
Genome annotation with 2820 genetic markers for budding yeast strains BY and RM and their progeny, 124 genotypes in total. Covers all the strains mentioned in brem2005_strains.csv and provides additional ones.

## Foss2007_protein_expression.txt:
A file from Zia Khan (June 2010) with protein abundance data. 574 columns. (Let's start numeration from 0.) Line 0: column names. Columns 0-16: additional data for a line. Column 1: protein name. Columns 17-464: values from 2 parents and 107 progeny, 10 replicates for each parent and 4 replicates for each progeny. Columns 17-444: progeny values. Columns 445-454: BY parent. Columns 455-464: RM parents. Columns 465-574: median values, for each progeny and for each parent. Columns 465-571: progeny median values. Column 572: BY median value. Column 573: RM median value.

## GDS1115.soft:
Genetic variation in gene expression among parents and progenies from brem2005_strains.csv
Expression profiling of parental strains and haploid progenies from a cross of strain BY4716 and the wild wine strain RM11-1a. Results used to find linkage between gene expression levels, which are treated as quantitative traits, and genetic markers.

## GDS1116.soft:
The same thing, but with dye-swap arrays.
