# loading required libraries
library(tidyverse)
library(reshape2)
library(Hmisc)

# calculate similarity
diff_expr <- file.path('results/diff_expr')
sim_dir <- file.path('results/similarity')
if (!dir.exists(sim_dir)) dir.create(sim_dir)

tibble(files = list.files(diff_expr, recursive = TRUE, full.names = TRUE)) %>%
    mutate(pert_type = str_split(files, '/|\\.tsv', simplify = TRUE)[, 3]) %>%
    group_split(pert_type) %>%
    map(function(trt) {
        trt_dir <- file.path(sim_dir, unique(trt$pert_type))
        if (!dir.exists(trt_dir)) dir.create(trt_dir)
        
        map(trt$files, function(x) {
            r <- read_tsv(x) %>%
                acast(ID~pert_iname, value.var = 'logFC') %>%
                rcorr()
            
            df <- left_join(melt(r$r, value.name = 'SCC'),
                            melt(r$P, value.name = 'PVAL')) %>%
                as_tibble() %>%
                na.omit() %>%
                mutate(cell_id = str_split(x, '/|\\.tsv', simplify = TRUE)[, 4],
                       pert_type = str_replace_all(unique(trt$pert_type), '\\.', ' ')) %>%
                mutate_if(is.numeric, round, digits = 4) %>%
                dplyr::select(pert_type, Var1, Var2, cell_id, everything())
            
            write_tsv(df, file.path(trt_dir, paste0(unique(df$cell_id), '.tsv')))
        })
    })
