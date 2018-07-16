Foss2007.txt :

a file from Zia Khan (June 2010) with protein abundance data. 574 columns. (Let's start numeration from 0.) Line 0: column names. Columns 0-16: additional data for a line. Column 1: protein name. Columns 17-464: values from 2 parents and 107 progeny, 10 replicates for each parent and 4 replicates for each progeny. Columns 17-444: progeny values. Columns 445-454: BY parent. Columns 455-464: RM parents. Columns 465-574: median values, for each progeny and for each parent. Columns 465-571: progeny median values. Column 572: BY median value. Column 573: RM median value.


foss.pic :

Foss 2007 dataset as ProteinAbundanceProfiles object. 1073 genes, each has 109 values. None values are allowed. Obtained from columns 465-574 of file Foss2007.txt (which was obtained from Zia Khan in June 2010). There are 107 progeny, and parents BY and RM.


fossBY.pic :

BY parent values from Foss 2007 dataset as ProteinAbundanceProfiles object. 1073 genes, each has 10 values. None values are allowed. Obtained from columns 445-454 of file Foss2007.txt (which was obtained from Zia Khan in June 2010). These are 10 replicates of parent BY.


fossRM.pic :

RM parent values from Foss 2007 dataset as ProteinAbundanceProfiles object. 1073 genes, each has 10 values. None values are allowed. Obtained from columns 455-464 of file Foss2007.txt (which was obtained from Zia Khan in June 2010). These are 10 replicates of parent RM.


brem.pic :

Brem 2005 expression dataset as GeneExpressionProfiles object. 5610 genes, each has 115 values. None values are allowed. Obtained from files GDS1115.soft and GDS1116.soft. These are two replicates of similar data: 6 replicates for BY parent, 12 replicates for RM parent, and 113 segregants. Take 113 segregants, and for each, average between two replicates (if only one value is available, take this value, if no value is available, take None). As two last coordinates, take values for BY and for RM parents as averages over all replicates from both files (12 and 24 values, respectively).


bremBY.pic :

BY parent values from Brem 2005 expression dataset as GeneExpressionProfiles object. 5610 genes, each has 12 values. None values are allowed. Obtained from files GDS1115.soft and GDS1116.soft. Take 6 replicates for BY parent from each file.


bremRM.pic :

RM parent values from Brem 2005 expression dataset as GeneExpressionProfiles object. 5610 genes, each has 24 values. None values are allowed. Obtained from files GDS1115.soft and GDS1116.soft. Take 12 replicates for RM parent from each file.


complexes.pic :

the list of 600 complexes as MacromolecularComplex objects: 380 complexes from go_protein_complex_slim.tab file from yeastgenome.org and 220 complexes from mips.txt (obtained from Jimin Song in May 2010).


geneids.pic :

a GeneIdentifiers object obtained from file BIOGRID-IDENTIFIERS-3.0.66.tab.txt after filtering using FilterBiogridIdentifiersFile.


ppi.pic :

PPI network as a Biograph object. Obtained from BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.0.65.tab2.txt by a number of filterings. Vertices and edges are annotated with foss and brem data.