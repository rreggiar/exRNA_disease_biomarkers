#!/public/groups/kimlab/.install_bin/anaconda3/envs/aale.analysis.env/bin/Rscript

library(tidyverse)

paths.in <- scan(file=file("stdin", "r"), what="character", n=2)

input.file <- file.path(paths.in[1])
output.file <- file.path(paths.in[2])

aale.locus.count <-
	read_csv(input.file) %>%
	filter(grepl('Alu', X1)) %>%
	filter_at(vars(contains('.')), any_vars(. >= 10))

write_csv(aale.locus.count, output.file)
