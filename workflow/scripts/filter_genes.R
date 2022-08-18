library(dplyr)
library(tidyr)
library(seqinr)
library(bioformatr)

input  <- snakemake@input
output <- snakemake@output
params <- snakemake@params

tsv_file    <- unlist(input["tsv"])
fasta_files <- unlist(input["fasta"])
pfam_files  <- unlist(input["pfam"])
output_dir  <- unlist(output["dir"])
output_tsv  <- unlist(output["tsv"])

min.species <- unlist(params["min_species"])

dir.create(output_dir)

fasta <- lapply(fasta_files, read.fasta, as.string = T, seqtype = "AA") %>%
    lapply(function(x) data.frame(protein_id = names(x), transl = unlist(x), annot = sapply(x, attr, "Annot"))) %>%
    setNames(basename(fasta_files)) %>%
    bind_rows(.id = "fname") %>%
    mutate(genome = gsub(".faa", "", fname)) %>%
    mutate(annot = gsub("^[^ ]+ |putative ", "", annot)) %>%
    mutate(annot = ifelse(grepl("hypothetical protein|unknown|EsV-1|membrane protein", annot), NA, annot))

col_names <- c("protein_id", "aln_start", "aln_end", "env_start", "env_end", "hmm_acc", "hmm_name", "type", "hmm_start", "hmm_end", "hmm_length", "bit_score", "E_value", "significance", "clan")
pfam_scan <- lapply(pfam_files, read.table, col.names = col_names) %>%
    bind_rows
helios <- filter(pfam_scan, hmm_name == "Heliorhodopsin") %>%
    pull(protein_id)

genes <- read.table(tsv_file, comment.char = "", sep = "\t", header = T, na.strings = "*", check.names = F) %>%
    mutate(group = 1:n()) %>%
    filter(`# Species` == Genes, Genes >= min.species) %>%
    gather(fname, protein_id, contains(".faa"), na.rm = T) %>%
    filter(! protein_id %in% helios) %>%
    mutate(genome = gsub(".faa", "", fname)) %>%
    select(-fname) %>%
    left_join(fasta, by = c("genome","protein_id")) %>%
    group_by(group) %>%
    group_walk(function(x, y) {
        fname <- sprintf("%s/%03d.faa", output_dir, y$group)
        mutate(x, sprintf(">%s %s %s\n%s", genome, protein_id, ifelse(is.na(annot), "", annot), transl)) %>%
            pull %>%
            cat(file = fname, sep = "\n")
    })
group_annot <- filter(genes, !is.na(annot)) %>%
    group_by(group, annot) %>%
    summarize(n = n(), annot = sprintf("%s (%d)", first(annot), n), .groups = "drop_last") %>%
    arrange(-n) %>%
    summarize(annot = paste(annot, collapse = "; "), .groups = "drop")
group_pfam <- left_join(pfam_scan, genes, by = "protein_id") %>%
    mutate(hmm = paste(hmm_acc, hmm_name, sep = " - ")) %>%
    group_by(group, hmm) %>%
    summarize(n = n(), pfam = sprintf("%s (%d)", first(hmm), n), .groups = "drop_last") %>%
    arrange(-n) %>%
    summarize(pfam = paste(pfam, collapse = "; "), .groups = "drop")

select(genes, group, `# Species`, Genes, `Alg.-Conn.`, genome, protein_id) %>%
    left_join(group_annot, by = "group") %>%
    left_join(group_pfam, by = "group") %>%
    spread(genome, protein_id) %>%
    rename(Orthogroup = group, Annotation = annot, Pfam = pfam) %>%
    write.table(file = output_tsv, row.names = F, sep = "\t", na = "", quote = F)
