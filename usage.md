
# D-BPP Pipeline: Step-by-Step Workflow

This document provides a detailed, step-by-step guide for running the D-BPP pipeline, from raw data preparation to the final inference of phylogenetic networks and ghost introgression.

---
## 1. Software Requirements
### a. Dsuite: https://github.com/millanek/Dsuite/tree/master
Citation:

Malinsky, M., Matschiner, M. and Svardal, H. (2021) Dsuite ‐ fast D‐statistics and related admixture evidence from VCF files. Molecular Ecology Resources 21, 584–595. doi: https://doi.org/10.1111/1755-0998.13265

```
git clone https://github.com/millanek/Dsuite.git
cd Dsuite
make
```
### b. BPP: https://github.com/bpp/bpp
Citation:

Flouri T., Jiao X., Rannala B., Yang Z. (2018) Species Tree Inference with BPP using Genomic Sequences and the Multispecies Coalescent. Molecular Biology and Evolution, 35(10):2585-2593. doi:10.1093/molbev/msy147
```
wget https://github.com/bpp/bpp/releases/download/v4.8.4/bpp-4.8.4-linux-x86_64.tar.gz
tar zxvf bpp-4.8.4-linux-x86_64.tar.gz
```
### c. snp-sites: http://sanger-pathogens.github.io/snp-sites/
Citation：
Page AJ, Taylor B, Delaney AJ, Soares J, Seemann T, Keane JA, Harris SR. 2016. SNP-sites: Rapid efficient extraction of SNPs from multi-fasta alignments. Microb Genom 2:e000056

Install Conda and install the bioconda channels.
```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install snp-sites
```

## 1. Data Preparation

### a. Sequence Data

- **Input:** Multi-locus sequence alignments, including both ingroup species and at least one designated outgroup.
- **Format:** Multi-FASTA. Ensure sample names are consistent throughout all analyses.

#### Example: Multi-FASTA Format

```fasta
>Sp1
ATGCTAG...
>Sp2
ATGCCAG...
>Sp3
...
>Outgroup
...

