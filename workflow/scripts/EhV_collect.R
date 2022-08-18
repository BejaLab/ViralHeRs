
library(dplyr)
library(tidyr)

input <- snakemake@input
output <- snakemake@output

blast_files <- unlist(input)
output_file <- unlist(output)

non_empty_files <- blast_files[file.size(blast_files) > 0]

lapply(non_empty_files, read.table, sep = "\t") %>%
    setNames(basename(non_empty_files)) %>%
    bind_rows(.id = "fname") %>%
    select(fname, qseqid = 2, sseqid = 3, pident = 4, stitle = 5) %>%
    distinct(fname, qseqid, .keep_all = T) %>%
    extract(stitle, into = "clade", regex = "clade=([^ ]+)") %>%
    mutate(virus = gsub("[.].+", "", fname)) %>%
    select(virus, qseqid, clade) %>%
    spread(clade, qseqid) %>%
    write.table(file = output_file, sep = "\t", row.names = F, quote = F)
