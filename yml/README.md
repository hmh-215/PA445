# Conda Environments for PA445

This directory contains conda environment YAML files specifying all environments and tools required by the workflows in this repository.

## Environment Inventory

| Environment File | Environment Name | Tools & Key Packages | Purpose |
|---|---|---|---|
| [`BPannotation.yml`](file:///E:/huong/Projects/github/PA445/yml/BPannotation.yml) | `BPannotation` | Bakta (1.12.0), Prokka (1.15.6), ABRicate (1.4.0), eggNOG-mapper (2.1.13) | Bacterial annotation & AMR detection |
| [`BPstructure.yml`](file:///E:/huong/Projects/github/PA445/yml/BPstructure.yml) | `BPstructure` | CheckM (1.2.5), CheckV (1.0.3) | Completeness & contamination QC |
| [`BPtyping.yml`](file:///E:/huong/Projects/github/PA445/yml/BPtyping.yml) | `BPtyping` | mlst (2.23.0), pasty | MLST and serotyping |
| [`assembly.yml`](file:///E:/huong/Projects/github/PA445/yml/assembly.yml) | `assembly` | SPAdes (4.0.0), QUAST (5.2.0), seqkit (2.13.0), Unicycler (0.5.1) | Assembly assessment & sequence manipulation |
| [`recombination.yml`](file:///E:/huong/Projects/github/PA445/yml/recombination.yml) | `recombination` | ISEScan (1.7.3), IntegronFinder (2.0.6), EMBOSS (6.6.0.0), PhiSpy (5.0.10) | IS elements, integrons, prophages |
| [`phylogeny.yml`](file:///E:/huong/Projects/github/PA445/yml/phylogeny.yml) | `phylogeny` | clinker (0.0.32), MAFFT (7.525), IQ-TREE (3.1.2) | Multiple sequence alignment, phylogeny, synteny |
| [`ncbi.yml`](file:///E:/huong/Projects/github/PA445/yml/ncbi.yml) | `ncbi` | BLAST (2.12.0), NCBI datasets CLI (18.26.0), AMRFinderPlus (4.2.7), Python helper scripts | NCBI downloads, AMR profiling, Python helpers |
| [`rgi_env.yml`](file:///E:/huong/Projects/github/PA445/yml/rgi_env.yml) | `rgi_env` | RGI (6.0.5) + CARD database | Comprehensive AMR gene screening |
| [`pangenome.yml`](file:///E:/huong/Projects/github/PA445/yml/pangenome.yml) | `pangenome` | Panaroo (1.6.0) | Core-genome extraction |
| [`ani.yml`](file:///E:/huong/Projects/github/PA445/yml/ani.yml) | `ani` | fastANI (1.34), aniclustermap (2.0.1) | Average Nucleotide Identity calculation |
| [`plasmid.yml`](file:///E:/huong/Projects/github/PA445/yml/plasmid.yml) | `plasmid` | Platon (1.7), PlasmidFinder (2.1.6), MOB-suite (3.1.9) | Plasmid identification and typing |
| [`beast_env.yml`](file:///E:/huong/Projects/github/PA445/yml/beast_env.yml) | `beast_env` | BEAST2, TreeAnnotator, TreeTime, snp-sites, Gubbins | Temporal phylogenetic dating |
| [`snippy_env.yml`](file:///E:/huong/Projects/github/PA445/yml/snippy_env.yml) | `snippy_env` | Snippy | Core SNP calling against reference |
| [`r_env.yml`](file:///E:/huong/Projects/github/PA445/yml/r_env.yml) | `r_env` | R, optparse, ggtree, treeio, ggplot2, tidyverse | Tree and heatmap plotting figures |

## Setup Instructions

To create an individual environment:
```bash
conda env create -f yml/<environment_file>.yml
```
