# D-BPP Pipeline: Step-by-Step Workflow

**Yang Y, Pang XX, Ding YM, Zhang BW, Bai WN, Zhang DY. 2025. Synergizing Bayesian and Heuristic Approaches: D-BPP Uncovers Ghost Introgression in Panthera and Thuja. bioRxiv 2025.06.27.662067. doi: https://doi.org/10.1101/2025.06.27.662067**


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

If you have multiple locus-specific FASTA files, you should concatenate them into a single multi-locus alignment.

```
perl concatenate_multi-locus.pl file_list indir output.fasta
```
**Parameters:**

- file_list: A text file containing the names of all .fa files to be concatenated, one filename per line.

- indir: Directory containing the input .fa files.

- output.fasta: Name of the merged output FASTA file.

Description:
concatenate_multi-locus.pl:This script reads a list of .fa files specified in file_list from the directory indir, extracts and concatenates sequences for each unique sequence ID across all files, and writes the merged results to output.fasta. The sequence IDs are normalized by removing special characters and extending/truncating to 14 characters for consistency.

Example:
Suppose you have a folder ./fa_dir with files sample1.fa, sample2.fa, sample3.fa, and a file_list.txt containing:
```
sample1.fa
sample2.fa
sample3.fa
...
```

Then convert the combined file to VCF format.

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
- **loci.bpp**         	          ## multi-locus file
- **imap.txt**                    ## individual-to-species mapping
- **bpp.ctl**                     ## parameter configuration file

### b. Prepare for loci.bpp File

First, a loci list file (e.g., loci.list) should be prepared, listing the names of all loci files (one per line). 
Example: loci.list
```
locus1
locus2
locus3
...

```

Store all loci files in the input folder with a .fa extension, then run the following command:

```
perl fasta2bpp.pl input loci.list > loci.bpp
```
The example of loci.bpp
```
6 100
loci1^A1  ATGCCGTAGCTAGTCA...GATC
loci1^A2  ATGCCGTAGCTAGTCA...GATC
loci1^B1  ATGTCGTAGCTAGTCA...GATC
loci1^B2  ATGTCGTAGCTAGTCA...GATC
loci1^C1  ATGCCGTAGATAGTCA...GATC
loci1^C2  ATGCCGTAGATAGTCA...GATC

6 95
loci2^A1  TCGTAGTCAAGGTAC...CTGA
loci2^A2  TCGTAGTCAAGGTAC...CTGA
loci2^B1  TCGTAGTCAAGGTGC...CTGA
loci2^B2  TCGTAGTCAAGGTGC...CTGA
loci2^C1  TCGTAGTCAAGGTAT...CTGA
loci2^C2  TCGTAGTCAAGGTAT...CTGA

```
**Explanation:**

- Each block represents one locus (gene/region).

- The first line of each block gives the number of samples and the sequence length for that locus (e.g., 6 100 means 6 samples, each with a 100 bp sequence).

- Each subsequent line provides a unique sample identifier, typically in the format locusID^sampleID (e.g., loci1^A1), followed by the aligned DNA sequence for that locus.

- The order and names of samples must be identical across all loci.

- Each locus block is separated by a blank line.

- All sequences within a locus must have equal length and be properly aligned.


### c. Prepare for imap.txt File

The imap file defines the correspondence between individual samples and species, with each line containing an individual and its assigned species separated by a tab character.

```
A1  A
A2  A
B1  B
B2  B
C1  C
C2  C
```
### d. Prepare for bpp.ctl File

The example of bpp.ctl

```
  seed =  -1
  seqfile = loci.bpp
  Imapfile = imap.txt
  outfile = bpp.out
  mcmcfile = bpp.mcmc
  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
* speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio
  speciestree = 0        * species tree fixed

* speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 6  A B C D E F
		    2 2 2 2 2 3
  (F, ((D)J[&phi=0.200000,tau-parent=yes],((K[&phi=0.800000,tau-parent=yes],(E,M[&phi=0.100000])N)X,(A,((B,N[&phi=0.200000])M,C)S)T)R)K)O;
  usedata = 1  * 0: no data (prior); 1:seq like
  nloci = 1000  * number of data sets in seqfile
  phase =   1 1 1 1 1 1
  cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

  thetaprior = 3 0.002 e  # invgamma(a, b) for theta
  tauprior = 3 0.01    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's
  phiprior = 1 1  # Beta(a, b) for root tau & Dirichlet(a) for other tau's

* heredity = 1 4 4
* locusrate = 1 15
  locusrate = 1 0 0 5 iid
  finetune =  0: 0.05 0.001 0.001 0.001 0.05 0.001 0.001 # finetune for GBtj, GBspr, theta, tau, mix
  Threads = 20 1 1
  print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
  burnin = 150000
  sampfreq = 2
  nsample = 1500000
```

Be sure to edit the input file parameters and the species and tree configurations according to your specific dataset. Additionally, update the designated introgression models as needed before running the analysis.

For more detailed settings, please refer to the official BPP manual.


### e. run BPP
```
bpp -cfile bpp.ctl
```
