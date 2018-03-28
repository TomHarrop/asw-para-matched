#!/usr/bin/env Rscript

library(data.table)

###########
# GLOBALS #
###########

data_dir <- "data/asw_para_matched"
log_file <- paste(data_dir, "combine.log", sep = "/")

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# find data files
key_files <- list.files(data_dir,
                        recursive = FALSE,
                        pattern = ".txt",
                        full.names = TRUE)
names(key_files) <- sub(".txt", "", basename(key_files))
print(key_files)

# combine datafiles
key_data <- rbindlist(lapply(key_files, fread),
                      idcol = "key")

# remove ophir samples
filtered_key_data <- key_data[!startsWith(sample, "O")]

# fix sample names
filtered_key_data[endsWith(sample, "*"),
                  sample_name := gsub("^R([[:digit:]]+).*", "Rpoa\\1", sample)]
filtered_key_data[endsWith(sample, "-H"),
                  sample_name := gsub("-H$", "", sample)]
filtered_key_data[startsWith(sample, "GBSNEG"),
                  sample_name := paste0(flowcell, lane, sample)]
filtered_key_data[is.na(sample_name), sample_name := sample]

# define populations
filtered_key_data[, population := gsub("[[:digit:]]+", "", sample_name)]
filtered_key_data[startsWith(sample, "GBSNEG"), population := "GBSNEG"]

# select columns
munged_data <- filtered_key_data[, .(
    key, flowcell, lane, barcode,
    agr_sample_name = sample,
    sample_name,
    population
)]

# write output
fwrite(munged_data, paste(data_dir, "combined_key_data.csv", sep = "/"))

# write session info
sessionInfo()
