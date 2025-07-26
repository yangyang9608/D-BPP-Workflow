
# D-BPP Pipeline: Step-by-Step Workflow

This document provides a detailed, step-by-step guide for running the D-BPP pipeline, from raw data preparation to the final inference of phylogenetic networks and ghost introgression.

---

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

