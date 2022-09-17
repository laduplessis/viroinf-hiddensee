# Datasets

### (VIROINF - Workshop on analysis of viral NGS data & evolutionary history and co-phylogenetics)

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
Dataset from Metzker et al ([PNAS 2002](http://www.pnas.org/content/99/22/14292)). Downloaded from [https://pubmed.ncbi.nlm.nih.gov/12388776/](https://pubmed.ncbi.nlm.nih.gov/12388776/) for both env (PopSet 24210023; https://www.ncbi.nlm.nih.gov/popset/24210023) and pol (PopSet 24209939; https://www.ncbi.nlm.nih.gov/popset/24209939) datasets.


### Env dataset 
Dataset of 132 gp120 sequences. Aligned in AliView using MUSCLE, trimmed ends and insertions in only in one sequence and renamed sequence headers. Resulting alignment is 820 bp long.


### Pol dataset
Dataset of 42 pol sequences. Aligned in AliView using MUSCLE, trimmed ends and renamed sequence headers. Resulting alignment is 689 bp long.


### Sequence headers
In format `<Genbank id>|<name>`.

- Victim sequence names start with V
- Patient sequence names start with P
- Control group (from the Lafayette area) sequence names start with LA

For technical replicates, sequences containing BCM in the name were sequenced at Baylor College of Medicine and sequences with MIC in the name were sequenced at the University of Michigan. 