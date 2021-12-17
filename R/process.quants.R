#!/user/bin/env Rscript

## Roman Reggiardo
## 2019.07.04

# in: [1] Directory containing salmon output subdirectories with naming (ctrl.1,
# kras.2.....)
#     [2] Sleuth LRT output file path
#
# out: Collapses all quant.sf files and calculates log2(Fold Change) values from TPM
# in each condition ['quant.data'], merges with the Sleuth LRT output to
# generate a comprehensive tx, p/q value, log2fc data set ['quant.data.lrt']


#####
library(tidyverse)
#####

# read stdin length of 2
paths.in <- scan(file=file("stdin", "r"), what="character", n=2)
# set first stdin argument to paths (will be iterated over)
paths <- paths.in[1]
# set second stdin argument to sleuth output file
sleuth.out <- read.csv(paths.in[2])
#sample.names <- tibble(sample = c('ctrl.1', 'ctrl.2', 'ctrl.3',
#                                  'kras.1', 'kras.2', 'kras.3'))
sample.names <- tibble(sample = c('intra.kras.1','intra.kras.2','intra.kras.3',
                                  'exo.kras.1', 'exo.kras.2', 'exo.kras.3'))
# extract subdir's from input path
files <- file.path(paths, sample.names$sample, 'te.locus.out.for.sleuth')
# print subdir names
print(files)
# set files list names to the names of the subdir's 
names(files) <- files
# [iteration]
files <- lapply(files, function(x){
                 # extract contents of subdirs by pasting together subdir
                 # and upstream path     
                      full.path <- list.files(path=x,
                                                  pattern='*sf', full.name = T)
                      # read in quant.sf files
                      x <- read.delim(full.path, sep = '\t')
                      # set new rep column to name of subdir (condition)
                      print(full.path)
                      x$rep <- strsplit(full.path,'/')[[1]][11]
                      print(strsplit(full.path,'/')[[1]][11])
                      # seperate rep column by the file name delimeter into
                      # condtion 'ctrl/kras' and rep '1,2,3'
                      x <- x %>% separate(rep, sep = '[.]', c('condition', NA, 'rep'))
                      # return the new dataframe
                      return(x)
}) # end [iteration]
# bind all of the data frames created above row-wise
quant.data <- bind_rows(files)
# [pipe]
head(quant.data)
quant.data <- quant.data %>% 
# group by tx name and condition(kras/ctrl)
group_by(Name, condition) %>% 
# take the tx,condition average TPM
summarize(avg = mean(TPM)) %>%
#mutate(avg = avg + 1) %>%
# take the condition entries and make them columns with the values set to TPM
spread(condition, avg)
# calculate lof2(Fold Change) in a new column
quant.data$log2fc <- log2((quant.data$exo + 1.00)/(quant.data$intra + 1.00))
# merge the sleuth output with the new quant data by the tx names
quant.data.lrt <- merge(quant.data, sleuth.out, by.x = 'Name', by.y = 'target_id')
# preview output file
head(quant.data.lrt)
# save quant.data.lrt as .csv in subdir's directory (may want to change this)
write.csv(quant.data.lrt, paste(paths, 'quant.data.log2fc.csv', sep = '/'))
