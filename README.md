# sarscov2-tools
Scripts for analysing SARS-CoV-2 nucleotide sequences

This repository contains scripts that are useful for calculating average genetic distance, wu-kabat variability coefficient and lineage frecuency per 2020 week, among others

### Python packages
- pandas
- numpy
- argparse
- collections
- glob
- re
- os


## Average genetic distance
[average_distance.py](https://github.com/mlarjim/sarscov2-tools/blob/main/average_distance.py/) determines average genetic distance and standard deviation of a symmetric matrix. This matrix must be generated with [snp-dists](https://github.com/tseemann/snp-dists) (pairwise SNP distance matrix from a FASTA sequence alignment). snp-dists output is average_distance.py input.
Usage:

```
average_distance.py snp-dists-matrix
```


## Coexisting lineages per week
[coexisting_linages.py](https://github.com/mlarjim/sarscov2-tools/blob/main/coexisting_linages.py/) outputs a table (.tsv) with the frecuency of every lineage per week.

```
coexisting_lineages.py metadata [-f] [-p]

metadata  tsv with (at least) 'Lineage' and 'date' columns
-f  outputs .tsv with the lineages' earliest date
-p  outputs .tsv of lineages per Andalusian province
```


## Mutations frecuency per week
[cov_cluster.py](https://github.com/mlarjim/sarscov2-tools/blob/main/cov_cluster.py/) outputs a table with the frecuency per week of a specific amino acid mutation of protein S.

```
cov_cluster.py metadata directory1 position

metadata  tsv with (at least) 'Strain' and 'date' columns
directory1  /complete/directory/to/*_analysis_report.csv
position  comma separated amino acids positions of protein S
```


## Mutations frecuency per province
[cov_provincias.py](https://github.com/mlarjim/sarscov2-tools/blob/main/cov_provincias.py/) outputs a table with the frecuency of a specific mutation of protein S per week per Andalusian province.

```
cov_provincias.py metadata directory1 position [-n]

metadata  tsv with (at least) 'Strain', 'date' and 'Province' columns
directory1  /complete/directory/to/*_analysis_report.csv
position  comma separated amino acid position of protein S
-n  outputs table with total number of samples per week per province
```


## Mutations per sample and frecuency of mutations
[frec_snp_gtc.py](https://github.com/mlarjim/sarscov2-tools/blob/main/frec_snp_gtc.py/) performs variant calling from .tsv files (output of [iVar](https://github.com/andersen-lab/ivar). Only analyse mutations with a frecuency >=0.75 of sequences with a coverage >=0.95.

```
frec_snp_gtc.py  directory1 [-t]
directory1  /complete/directory/to/*_analysis_report.csv
-t  outputs .tsv file of nucleotide substitution type frecuency (for example, C>T frecuency)
```


## Synonymous and non-synonymous mutations
[tabla_Zekri_aa.py](https://github.com/mlarjim/sarscov2-tools/blob/main/tabla_Zekri_aa.py/) analyse mutations per gene of SARS-CoV-2.

```
tabla_Zekri_aa.py  directory1 [-a]
directory1  /complete/directory/to/*_analysis_report.csv
-a  adds a columns of amino acids alterations to .tsv file 
```


## Wu-Kabat variability coefficient
[wu_kabat_gtc.py](https://github.com/mlarjim/sarscov2-tools/blob/main/wu_kabat_gtc.py/) calculates Wu-Kabat variability coefficient and outputs a .tsv file for every protein.

```
wu_kabat_gtc.py  directory1 [-m] [-p]
directory1  /complete/directory/to/*_analysis_report.csv
-m  outputs a single .tsv file with all amino acid mutations 
-p  outputs a .tsv file for every protein with variability coefficient even if there is no variability in that amino acid.
```

## Bibliography
Wu, T. T. & Kabat, E. A. An analysis of the sequences of the variable regions of Bence Jones proteins and myeloma light chains and their implications for antibody complementarity. J. Exp. Med. 132, 211–250 (1970).

Zekri, A. N. et al. Genomic characterization of SARS-CoV-2 in Egypt. J. Adv. Res. (2020).



#### Owner note
Some scripts are aimed to Andalusian (from Andalusia, Spain) sequences because they were created as part of my Thesis proyect: *Genomic and population diversity analysis of SARS-CoV-2 during first epidemic wave in Andalusia*

🍃

LinkedIn [🔗](www.linkedin.com/in/maria-lara-jimenez)

