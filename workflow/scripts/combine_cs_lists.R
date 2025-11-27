#!/usr/bin/Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

#----------#
# loading params
sets_path <- snakemake@input
file_path_res <- snakemake@output[["cslist"]]
file_path_rep <- snakemake@output[["report"]]

#--------------#
# results path
exmpl_path <- as.character(sets_path[1])
pcmd <- dirname(exmpl_path)

#--------------#
# scan results for all proteins sequence
res_files <- list.files(
  pattern = paste0("seq.(\\d+).(\\d+)_(\\d+)_(\\d+)_(\\d+).cslist"),
  path = pcmd,      # the path where the independents snps file live
  recursive = FALSE, # to show the files in subdirectories or subfolders
  full.names = TRUE # to show full path
)

#--------------#
# extract seqid from input sentinel files
input_seqid <- map_dfr(
  sets_path, function(path) {
      base_path = dirname(path)
      file_name = basename(path)
      seqid = gsub("_.*$", "", file_name) # extract the protein sequence id and remove file format
  data.frame(base_path, file_name, seqid)
  }
)

#--------------#
# extract seqid from COJO outputs
seq_list_tbl <- tibble(res_files) %>% mutate(seqid = str_extract(res_files, "seq.\\d+.\\d+"))

# select input seqids from the existing results
res_files_input <- res_files[seq_list_tbl$seqid %in% input_seqid$seqid]

# report files based on input results file names
rep_files_input <- res_files_input %>%
  str_replace("cs_list", "cs_report") %>% # rename folder, then rename file format
  str_replace(".cslist", ".report")

# combine results for the input seqids
res_combined <- map_dfr(res_files_input, fread)
rep_combined <- map_dfr(rep_files_input, fread)

#--------------#
# save the joint results
write.table(res_combined, file = file_path_res, sep = "\t", quote = F, row.names = F)
write.table(rep_combined, file = file_path_rep, sep = "\t", quote = F, row.names = F)
