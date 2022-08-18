Viral heliorhodopsins pipeline
==============================

This repository contains the [Snakemake](https://snakemake.readthedocs.io/) workflow used for the bioinformatic analyses for the paper Hososhima et al (2022) "[Proton-transporting heliorhodopsins from marine giant viruses](https://doi.org/10.7554/eLife.78416)" eLife.

All of the dependencies are taken care of with conda, so it is recommended to run `snakemake` with `--use-conda`.

This repository is organized as follows:

* `analysis` -- contains the intermediate files
* `annotations` -- manually curated data
* `databases` -- includes Pfam databases and algal protein sequences. To run the workflow from scratch, the soft link `databases/Pfam` should point to the Pfam database folder and soft links in `databases/algae` should point the corresponding fasta files.
* `output` -- final output files
* `proteins` -- curated sequence data for algal and viral heliorhodopsins
* `viruses` -- GenBank files with the viral genomes
* `workflow` -- workflow files, including:
    * `envs` -- conda environment files
    * `Snakefile` -- the snakemake file
    * `scripts` -- folder with scripts

The output files are as follows:

* `cat_phylogeny.pdf` -- concatenation phylogeny of the viruses
* `EhVHeRs.tsv` -- distribution of heliorhodopsin genes among EhVs
* `HeR_tree.jtree` -- phylogenetic tree of viral and algal heliorhodopsins in `.jtree` format
* `HeR_tree.pdf` -- image version of the same tree
* `miniset_chronos.pdf` -- small HeR tree with alignment of critical positions
* `orthogroups.tsv` -- EhV orthogroups
