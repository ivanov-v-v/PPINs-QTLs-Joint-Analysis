# Mutual analysis of interaction networks and quantitative trait loci for budding yeast

---
## Table of contents
* [Description](#description)
* [Motivation](#motivation)
* [How to use this?](#how-to-use-this)

---
## Description

![Saccharomyces cerevisiae, the budding yeast](http://ppdictionary.com/mycology/budding_yeast.jpg)
This is my first project in systems biology. I started it during the second year of my undergraduate studies in Moscow 
Institute of Physics and Technology ([MIPT](https://mipt.ru/english/), in Russia better known as "Физтех").  

This research was carried out under supervision of [Yuri Pritykin](https://scholar.google.com/citations?hl=en&user=Arx56RkJBrYC&view_op=list_works&sortby=pubdate).
This GitHub is aimed to ensure the reproducibility of our results.

[↥ back to top](#table-of-contents)

---
## Motivation

* **Why QTLs?** Quantitative trait loci or QTLs are a hot topic in contemporary biological research. The link between 
genotype and phenotype, they attracted increasing attention in recent years due to developments of high-throughouput 
technologies like RNA-seq. This research deals with eQTLs — with mRNA expression as quantitative trait — and pQTLs — 
same for protein expression. Genetic architectures of eQTLs and pQTLs are very different, despite the obvious 
relationship between transcription and translation of the genes. There is an ongoing debate on this matter in science.
We decided to make a contribution by comparing genetic architectures underlying eQTL and pQTL sets linked to functional
modules — sets of interacting genes participating in a particular cell process. Our approach is the possible solution 
that reduces eQTL/pQTL discrepancy via protein-protein interaction networks.          
* **Why yeast?** Despite our hypothesis isn't specific to this species, we concentrated on *Saccharomyces cerevisiae* 
or *baker's yeast*. It's a popular model organism, sequenced to unprecedented depth. Most of its genes are 
annotated and their interactions are thoroughly studied. Yeast biological databases are well-developed, carefully 
curated and easy to work with, which results in increasing quality of computational experiments.  
* **Why does this matter?** Some complex diseases can be tracked down to genetic level — either mRNA or protein levels
of particular genes are wrong. Understanding the relationships between these two is necessary for the accurate treatment. 
eQTLs have been studied for the long time already, experiment designs are standardized and well-known. 
pQTLs are a different matter, because protein expression measurements are still costly and carry a lot
of noise. Introduction of additional structure — the interactome — helps us to make the most out of limited data. 
It may be, that to bridge the gap between eQTLs and pQTLs we should look at them on a bigger scale.             
* **Why us?** We have an unusually rich protein expression dataset. It drastically increases the statistical 
significance of our observations.

[↥ back to top](#table-of-contents)

---
## How to use this?

This is my first medium-scale project as well as the first research project, so please just bear with me and be patient.
Major refactoring is scheduled for the following month, it will get better soon.

This section will be subdivided into:
1. Instructions on how to pull this repository and how to configure the right conda environment. Maybe, this repository
will get wrapped in a Docker container in the future.
2. Description of directory structure.
3. Driver scripts that recompute certain parts of the project together with the order in which they
are supposed to be ran.
4. Link to documentation (Doxygen?)

[↥ back to top](#table-of-contents)
---