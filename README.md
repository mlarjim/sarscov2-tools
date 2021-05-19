# sarscov2-tools
Scripts for analysing SARS-CoV-2 nucleotide sequences

This repository contains scripts that are useful for calculating average genetic distance, wu-kabat variability coefficient and lineage frecuency per 2020 week, among others

### Python packages
- pandas
- numpy
- argparse
- collections

## Average genetic distance
[average_distance.py](https://github.com/mlarjim/sarscov2-tools/blob/main/average_distance.py/) determines average genetic distance and standard deviation of a simmetric matrix. This matrix must be generated with snp-dists (pairwise SNP distance matrix from a FASTA sequence alignment). [snp-dists](https://github.com/tseemann/snp-dists) output is average_distance.py input.
Usage:

```
average_distance.py snp-dists-matrix
```

## Coexisting lineages per week

