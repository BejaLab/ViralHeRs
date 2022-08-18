library(dplyr)
library(tidyr)
library(phytools)
library(ggtree)
library(ape)
library(seqinr)
library(treeio)
library(ggplot2)
library(ggnewscale)

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

tree_file <- unlist(input["tree"])
fasta_file <- unlist(input["fasta"])
genera_file <- unlist(input["genera"])
clades_file <- unlist(input["clades"])
positions_file <- unlist(input["positions"])
colorscheme_file <- unlist(input["colorscheme"])
output_pdf_file <- unlist(output["pdf"])
output_jtree_file <- unlist(output["jtree"])

res_colors <- read.table(colorscheme_file, comment.char = "", header = T, fill = T, sep = "\t") %>%
    separate_rows(Residues, sep = "") %>%
    filter(Residues != "") %>%
    with(setNames(Color, Residues))

positions <- read.table(positions_file, header = T)
ref.seq <- names(positions)

clades <- read.table(clades_file, header = T, comment.char = "") %>%
    arrange(Clade)
genera <- read.table(genera_file, header = T)

fasta <- treeio::read.fasta(fasta_file) %>%
    lapply(as.character) %>%
    lapply(toupper) %>%
    bind_cols
sequences <- seqinr::read.fasta(fasta_file, seqtype = "AA", as.string = T) %>%
    {data.frame(label = names(.), sequence = as.character(.))}
aln <- setNames(fasta, sub(" .+", "", names(fasta))) %>%
    filter(!!as.symbol(ref.seq) != "-") %>%
    filter(1:n() %in% positions[[1]]) %>%
    t %>% as.data.frame %>%
    setNames(positions[[1]])

taxa <- data.frame(header = names(fasta)) %>%
    extract(header, into = c("label", "Species"), remove = F, regex = "^(.+?) .+\\[(.+)\\]$") %>%
    extract(Species, into = "Genus", remove = F, regex = "^([^ ]+)") %>%
    left_join(genera, by = "Genus") %>%
    left_join(clades, by = "Clade")

tree <- read.tree(tree_file) %>%
    midpoint.root %>%
    as.treedata %>%
    as_tibble %>%
    left_join(taxa, by = "label") %>%
    left_join(sequences, by = "label") %>%
    mutate(full_label = ifelse(grepl("HeR clade", header), sprintf("%s [%s]", label, Species), Species)) %>%
    as.treedata
write.jtree(tree, file = output_jtree_file)

colors <- with(clades, setNames(Color, Clade))

p <- ggtree(tree) +
    geom_point2(aes(subset = !isTip & as.numeric(label) >= 90, size = as.numeric(label) >= 95, x = branch)) +
    geom_tiplab(aes(color = Clade, label = full_label), size = 2, align = T) +
    geom_treescale(width = 0.5) +
    scale_size_manual(values = c(1,2)) +
    scale_color_manual(values = colors, na.value = "black") +
    hexpand(0.2)
p <- gheatmap(p + new_scale_fill(), aln, offset = 3, width = 0.2) +
    scale_fill_manual(values = res_colors)

ggsave(output_pdf_file, p, width = 10, height = 10)
