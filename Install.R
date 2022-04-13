R -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org", dependencies = TRUE)'

R -e 'BiocManager::install(c("rtracklayer", "eisaR", "tximeta", "BSgenome"), force = TRUE, update = TRUE, ask = FALSE)'
