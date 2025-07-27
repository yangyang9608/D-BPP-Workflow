
# D-BPP Pipeline: Step-by-Step Workflow

This document provides a detailed, step-by-step guide for running the D-BPP pipeline, from raw data preparation to the final inference of phylogenetic networks and ghost introgression.

---
## 1. Software Requirements
### a. Dsuite: https://github.com/millanek/Dsuite/tree/master
Citation: Malinsky, M., Matschiner, M. and Svardal, H. (2021) Dsuite ‐ fast D‐statistics and related admixture evidence from VCF files. Molecular Ecology Resources 21, 584–595. doi: https://doi.org/10.1111/1755-0998.13265

```
git clone https://github.com/millanek/Dsuite.git
cd Dsuite
make
```
### b. BPP: https://github.com/bpp/bpp
Citation: Flouri T., Jiao X., Rannala B., Yang Z. (2018) Species Tree Inference with BPP using Genomic Sequences and the Multispecies Coalescent. Molecular Biology and Evolution, 35(10):2585-2593. doi:10.1093/molbev/msy147
```
wget https://github.com/bpp/bpp/releases/download/v4.8.4/bpp-4.8.4-linux-x86_64.tar.gz
tar zxvf bpp-4.8.4-linux-x86_64.tar.gz
```
### c. snp-sites: http://sanger-pathogens.github.io/snp-sites/
Citation: Page AJ, Taylor B, Delaney AJ, Soares J, Seemann T, Keane JA, Harris SR. 2016. SNP-sites: Rapid efficient extraction of SNPs from multi-fasta alignments. Microb Genom 2:e000056

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

```
>Sp1
ATGCTAG...
>Sp2
ATGCCAG...
>Sp3
AATCCTG...
>Sp4
AAGATTC...
……
>Outgroup
...
```
### b. Species Tree Topologies
Input: Candidate species tree topologies in Newick format (e.g., generated from previous studies or other analyses).
Format example:
```
(((Sp1,Sp2),Sp3),Outgroup);
(((Sp1,Sp3),Sp1),Outgroup);
```

## 2. D-statistic (ABBA-BABA test)

Use Dsuite or an equivalent tool to identify potential introgression events.

If you have multiple locus-specific FASTA files, you should concatenate them into a single multi-locus alignment and then convert the combined file to VCF format.

```

snp-sites -v -o output.fasta input.vcf
```

Once the multi-locus FASTA files have been concatenated and converted to VCF format, you can use Dsuite to compute D-statistics for introgression analysis.

```
Dsuite Dtrios input.vcf imap --tree=TREE_FILE.nwk -o outputfilename
```
Output: outputfilename_tree.txt with D-statistic, Z-scores, p-values, and site pattern number of ABBA, BABA, BBAA.

Note: Only retain significant D-statistic signals (e.g., p < 0.01) for downstream validation.

## 3. Calculate the 𝐷𝑝
For each significant D-statistic result, compute the 𝐷𝑝 (Hamlin et al. 2020) using the specified formula to further quantify the extent of introgression.

```
D𝑝=(ABBA-BABA)/(BBAA+ABBA+BABA)

```
Then, rank the significant events by their 𝐷𝑝 values to determine the order in which they will be incorporated into the analysis.

## 4. Prepare the input files required for BPP analysis

### a. Required Input Files
i)   loci.bpp   ## multi-locus file
ii)  bpp.ctl    ## parameter configuration file
iii) imap.txt   ## individual-to-species mapping

### b. Prepare for loci.bpp File

If all loci files are stored in the input folder and each file has a .fa extension, you can use the following command.

In addition, a loci list file (e.g., loci.list) should be prepared, listing the names of all loci files (one per line). 
Example: loci.list
```
locus1
locus2
locus3
...

```


```
perl fasta2bpp.pl input loci.list > loci.bpp
```






