# Datasets

## EBOV dataset
Dataset of 15 Zaïre Ebola virus (EBOV) genomes taken from [http://beast.community/ebov_local_clocks.html](http://beast.community/ebov_local_clocks.html), where the original references for the genomes can be found.

Genomes were extracted from the BEAST 1.10 XML files using BEASTGen and then split into coding (cds) and non-coding (ig) regions using the gene annotations:

- NP:		470-2689
- VP35:	3129-4151
- VP40:	4479-5459
- GP:		6039-6923,6923-8068
- VP30:	8509-9375
- VP24:	10345-11100
- L:		11581-18219 

Instead of duplicating the frameshift position in GP a gap was inserted, to maintain the reading frame. 

### CDS 
Coding regions of the 15 EBOV genomes (14517 bp)

### IG
Non-coding regions of the 15 EBOV genomes (4449 bp)


## Dengue-4 dataset
Dataset of 17 pre-aligned DENV-4 envelope protein sequences (1485 bp). Originally from Lanciotti, Gubler & Trent ([Journal of Gen Vir 1997](https://www.microbiologyresearch.org/content/journal/jgv/10.1099/0022-1317-78-9-2279)) and featured in Jombart et al ([Mol Ecol Res 2017](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12676a)) and Drummond & Rambaut ([BMC Evol Biol 2007](http://dx.doi.org/10.1186/1471-2148-7-214)).

Dataset is included as an example dataset with [BEAST 2](https://beast2.org). Numbers after names are the collection years. Available in NEXUS and fasta formats. Sequnces divide into two clades:

- Clade I:
	- D4Philip64
	- D4Philip84
	- D4Philip56
	- D4Thai84
	- D4Thai78
	- D4SLanka78
	- D4Thai63
- Clade II:
	- D4Elsal94
	- D4PRico86
	- D4Tahiti85
	- D4Mexico84
	- D4Brazi82
	- D4ElSal83
	- D4NewCal81
	- D4Tahiti79
	- D4Indon77
	- D4Indon76



## HIV-1 criminal case dataset
Dataset from Metzker et al ([PNAS 2002](http://www.pnas.org/content/99/22/14292)). Downloaded from [https://pubmed.ncbi.nlm.nih.gov/12388776/](https://pubmed.ncbi.nlm.nih.gov/12388776/) for both env ([PopSet 24210023](https://www.ncbi.nlm.nih.gov/popset/24210023)) and pol ([PopSet 24209939](https://www.ncbi.nlm.nih.gov/popset/24209939)) datasets.


### Env dataset 
Dataset of 132 gp120 gene sequences. Aligned in AliView using MUSCLE, trimmed ends and insertions in only in one sequence and renamed sequence headers. Resulting alignment is 820 bp long.


### Pol dataset 
Dataset of 42 pol gene sequences. Aligned in AliView using MUSCLE, trimmed ends and renamed sequence headers. Resulting alignment is 689 bp long.


### Sequence headers
In format `<Genbank id>|<name>`.

- Victim sequence names start with V
- Patient sequence names start with P
- Control group (from the Lafayette area) sequence names start with LA

For technical replicates, sequences containing BCM in the name were sequenced at Baylor College of Medicine and sequences with MIC in the name were sequenced at the University of Michigan. 


## H3N2 dataset
Dataset of 86 aligned Influenza/A hemagglutinin (HA) gene sequences (1698 bp) from New York sampled over multiple flu seasons. The dataset is a subsample of the dataset used in Rambaut et al ([Science 2008](http://www.nature.com/nature/journal/v453/n7195/full/nature06945.html)).

## Alpha UK dataset
Two alignments of SARS-CoV-2 genomes from the UK. All genomes and associated 
metadata are available from https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/.

### Alpha
Subsample of all Alpha (B.1.1.7) genomes from the UK that were sequenced between 
1 August 2020 and 1 February 2021 (the oldest Alpha genome only dates from 20 September 2020). 
We sampled 7 genomes at random from all genomes sequenced in every epiweek during this 
period. For epiweeks where fewer than 7 genomes are available we used all 
available genomes. This resulted in 154 genomes. 

### Background
subsampling all non-Alpha genomes from the UK that were sequenced between 
1 August 2020 and 1 February 2021 (the oldest Alpha genome only dates from 20 September 2020). 
We sampled 6 genomes at random from all genomes sequenced in every epiweek during this 
period. For epiweeks where fewer than 6 genomes are available we used all 
available genomes. This resulted in 160 genomes. 
