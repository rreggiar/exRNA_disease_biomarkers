#!/user/bin/env Rscript
library(wasabi)
library(sleuth)
library(rhdf5)

### make this a for loop for all the metadata files in the directory ###
#### must be run in cmd line as 'echo 'path' 'out path' 'num samples' | Rscript sleuthRun.R
paths <- scan(file=file("stdin", 'r'), what='character', n=1)
path <- paths[1]
out <- '/public/groups/kimlab/aale.kras/data/bulk.rna.seq/ha1em/output'
sample_num <- 6

s2c <- data.frame('sample' = c('ctrl.1','ctrl.2','ctrl.3',
                           'kras.1', 'kras.2', 'kras.3'),
              'condition' = c('ctrl', 'ctrl', 'ctrl',
                              'kras', 'kras', 'kras'))
                        
#dir.create(out)
files <- file.path(path, s2c$sample, 
                   'te.locus.out.for.sleuth')
#print(files)
#sfdirs <- file.path(path, files)[1:sample_num]
prepare_fish_for_sleuth(files)

#sample_id <- files

#kal_dirs <- file.path(path, sample_id)


s2c <- dplyr::mutate(s2c, path = files[1:sample_num])
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~condition, "full")
so <- sleuth_fit(so, ~1, "reduced")
so <- sleuth_lrt(so, "reduced", "full")

sleuth_table <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 1.05)
write.csv(sleuth_significant, paste(out,
        '/locus.only.sleuth_lrt_out_qval_below_0.05_HIGH_CONFIDENCE', sep=''))
quit()
## pulled from https://groups.google.com/forum/#!topic/kallisto-sleuth-users/t7i6gldIMDc seems to address bug in sleuth
cond <- factor(s2c$condition)
cond <- relevel(cond,ref='ctrl')
md <- model.matrix(~cond,s2c)
colnames(md)[2] <- 'conditionkras'
so <- sleuth_prep(s2c, md, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, "full")
models(so)

soWT <- sleuth_wt(so, which_beta='conditionkras', which_model = 'full')

soWT_results <- sleuth_results(soWT, 'conditionkras', 'wt')
soWT_significant <- dplyr::filter(soWT_results, qval <= 0.05)
write.csv(soWT_significant, paste(out,
        '/combined.ref.sleuth_wt_out_qval_below_0.05_BE_CAREFUL', sep=''))
