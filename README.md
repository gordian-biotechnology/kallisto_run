# kallisto_run
## kallisto bustools pipeline with scvi and scanpy
Command line tool for single cell analysis. Builds kallisto indexes run kallisto bustools (Single cell qunat) and makes matrix and h5 files for scanpy/scvi analysis.

# Installation

```
git clone https://github.com/gordian-biotechnology/kallisto_run.git
cd kallisto_run
pip install .
```
## Requirements
1. [kallisto](https://pachterlab.github.io/kallisto/download)
2. [bustools](https://github.com/BUStools/bustools)
3. [pyTorch](https://pytorch.org/get-started/locally/)
4. [scvi](https://github.com/YosefLab/scVI)
5. [scanpy](https://icb-scanpy.readthedocs-hosted.com/en/latest/index.html)

# Usage
Help with input can be accessed from the command line:

```
> kallisto_run -h

[-n NAME] [-rr] [-v VERSION]

optional arguments:
-h, --help            show this help message and exit
-i INPUT, --input INPUT
    Path to directory with fastqs (default: None)
-o OUT, --out OUT     Path to output directory. (default: None)
-dr, --dry_run        Don't run just print out commands. (default: False)
-ref REF_DIR, --ref_dir REF_DIR
    Path to directory with references, should contain
    cdna.all.fa files. If not found will download by
    species, if .idx isn't also found it will build it.
    (default: None)
-s SPECIES, --species SPECIES
    Species of run, default is mouse. (default: mouse)
-n NAME, --name NAME  Name of to append to files. (default: )
-rr, --rerun
-v VERSION, --version VERSION
    Version of technology to provide to kallisto bus.
    (default: 10xv3)
```

Minimal example command:

```
kallisto_run -i /path/to/fastq_files -ref /path/to/human_kallisto_ref -s human -v 10xv3
```
If no reference directory is found or kallisto idx and gtf file are not found, they will be downloaded (ensembl) and the index will be built.
