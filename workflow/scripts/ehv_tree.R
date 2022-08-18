library(dplyr)
library(tidyr)
library(phytools)
library(ggtree)
library(ape)
library(treeio)

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

tree_file <- unlist(input["tree"])
fasta_files <- unlist(input["fasta"])
output_file <- unlist(output)

orthogroups <- lapply(fasta_files, read.fasta) %>%
    lapply(names) %>%
    unlist %>%
    sub(" .+", "", .) %>%
    table %>%
    data.frame %>%
    setNames(c("label", "genes"))

tree <- read.tree(tree_file) %>%
    midpoint.root %>%
    as.treedata %>%
    as_tibble %>%
    left_join(orthogroups, by = "label") %>%
    mutate(full_label = sprintf("%s (%s)", label, genes)) %>%
    as.treedata

p <- ggtree(tree) +
    geom_tiplab(aes(label = full_label), size = 5, align = T) +
    geom_point2(aes(subset = !isTip & as.numeric(label) > 95, x = branch)) +
    geom_treescale(width = 0.5) +
    hexpand(0.2)

ggsave(output_file, p, width = 5, height = 4)
