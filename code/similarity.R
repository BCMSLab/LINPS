# loading required libraries
library(tidyverse)
library(reshape2)
library(Hmisc)

# calculate similarity
perturbs <- c("trt_cp", "trt_sh", "trt_oe",  "trt_oe.mut")

map(perturbs,
    function(x) {
    d <- paste0('results/', x, '/diff_expr')
    cell <- str_split(list.files(d),
                      '\\.',
                      simplify = TRUE)[, 1]
    
    map(cell, function(c) {
        r <- read_tsv(file.path(d, paste0(c, '.tsv'))) %>%
            acast(ID ~ pert_iname, value.var = 'logFC') %>%
            rcorr()
            # rcorr(type = 'spearman') # TODO: choose a better measure of similarity
        
        sim <- left_join(melt(r$r, value.name = 'SCC'),
                  melt(r$P, value.name = 'PVAL')) %>%
            as_tibble() %>%
            na.omit() %>%
            mutate(cell_id = c)
        
        sim_dir <- paste0('results/', x,'/similarity')
        dir.create(sim_dir)
        
        write_tsv(sim,
                  file.path(sim_dir, paste0(c, '.tsv')))
        
    })
})
