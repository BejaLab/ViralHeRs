library(dplyr)
library(tidyr)
library(phytools)
library(ggtree)
library(treeio)
library(ape)
library(ggplot2)
library(ggnewscale)

input <- snakemake@input
output <- snakemake@output

tree_file <- unlist(input["tree"])
fasta_file <- unlist(input["fasta"])
outgroups_file <- unlist(input["outgroups"])

chronos_file <- unlist(output["chronos"])
phylo_file <- unlist(output["phylo"])
positions_file <- unlist(input["positions"])
colorscheme_file <- unlist(input["colorscheme"])

res_colors <- read.table(colorscheme_file, comment.char = "", header = T, fill = T, sep = "\t") %>%
    separate_rows(Residues, sep = "") %>%
    filter(Residues != "") %>%
    with(setNames(Color, Residues))

positions <- read.table(positions_file, header = T)
ref.seq <- names(positions)

outgroups <- read.fasta(outgroups_file) %>%
    names %>%
    sub(" .+", "", .)

fasta <- read.fasta(fasta_file) %>%
    lapply(as.character) %>%
    lapply(toupper) %>%
    bind_cols

deflines <- data.frame(header = names(fasta)) %>%
    extract(header, into = c("label", "skip", "Species"), remove = F, regex = "^([^ ]+)( .*\\[(.+)\\])?$") %>%
    extract(header, into = "alias", remove = F, regex = "=([^ ]+)") %>%
    select(-skip)

set.root.label <- function(my.tree) {
    root.node <- rootnode(my.tree)
    as.treedata(my.tree) %>%
        as_tibble %>%
        group_by(parent) %>%
        mutate(label = ifelse(parent == root.node & label == "", last(label), label)) %>%
        `class<-`(c("tbl_tree","tbl_df","data.frame")) %>%
        as.treedata
}

phylo.tree <- read.tree(tree_file) %>%
    reroot(getMRCA(., outgroups), position = 0.5, edgelabel = T)
chronos.tree <- chronos(phylo.tree, lambda = 3) %>%
    `class<-`("phylo")

tree <- set.root.label(chronos.tree) %>%
    as_tibble %>%
    left_join(deflines, by = "label") %>%
    mutate(full_label = sprintf("%s [%s]", ifelse(is.na(alias), label, alias), Species)) %>%
    as.treedata

aln <- setNames(fasta, sub(" .+", "", names(fasta))) %>%
    filter(!!as.symbol(ref.seq) != "-") %>%
    filter(1:n() %in% positions[[1]]) %>%
    t %>%
    `colnames<-`(positions[[1]]) %>%
    data.frame(id = rownames(.)) %>%
    gather(pos, res, -id) %>%
    mutate(pos = sub("X", "", pos) %>% as.numeric %>% as.factor)
p <- ggtree(tree) +
    geom_point2(aes(subset = !isTip & as.numeric(label) >= 90, size = as.numeric(label) >= 95, x = branch)) +
    geom_tiplab(aes(label = full_label), size = 4, align = T) +
    geom_treescale(width = 0.5) +
    scale_size_manual(values = c(1,2)) +
    # scale_color_manual(values = colors, na.value = "black") +
    hexpand(0.2)

p <- facet_plot(p + new_scale_fill(), panel = "aln", data = aln, geom = geom_label, aes(x = pos, label = res, fill = res), color = "white", fontface = "bold", family = "mono") +
    scale_x_discrete() +
    scale_fill_manual(values = res_colors, guide = "none") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 1)
    )

p.phylo <- ggtree(phylo.tree) +
    geom_tiplab(aes(label = label), size = 4, align = T) +
    geom_treescale(width = 0.5) +
    hexpand(0.2)

ggsave(chronos_file, p, width = 4, height = 4)
ggsave(phylo_file, p.phylo, width = 4, height = 8)
