# loading required libraries
library(tidyverse)
library(reshape2)
library(NPA)
library(NPAModels)

# load models
# TODO: apply calculation to all models
models <- load_models('Hs')
# devtools::install('data/NPAModels/')
# preprocessNetworks()
# models <- load_models('Hs')

# models <- map(list.files('data/NPAModels/data', pattern = '*.rda', full.names = TRUE),
#               function(x) {
#     env <- new.env()
#     load(x, envir = env)
#     env[[ls(env)]]
# })

perturbs <- c("trt_cp", "trt_sh", "trt_oe",  "trt_oe.mut")

map(perturbs, function(x) {
    d <- paste0('results/', x, '/diff_expr')
    cell <- str_split(list.files(d),
                       '\\.',
                       simplify = TRUE)[, 1]

    map(cell, function(c) {
        print(c)
        diff_expr <- read_tsv(file.path(d, paste0(c, '.tsv')))
        groups <- diff_expr %>%
            dplyr::select(nodeLabel = ID,
                   t,
                   foldChange = logFC) %>%
            as.data.frame() %>%
            with(split(., diff_expr$pert_iname)) %>%
            map(function(x) {
                filter(x, !duplicated(x$nodeLabel))
            })
        
        all(map(groups, ~all(sum(duplicated(toupper(.x))))) == 0)
        
        # npa_list
        npa_list_dir <- paste0('results/', x,'/npa_lists')
        dir.create(npa_list_dir)
        
        npa_list <- compute_npa_list(
            groups,
            models,
            verbose = TRUE)
        
        write_rds(npa_list, file.path(npa_list_dir, paste0(c, '.rds')))
        
        # bif
        bif_dir <- paste0('results/', x,'/bif')
        dir.create(bif_dir)
        
        bif <- get_bif(npa_list)
        
        write_rds(bif, file.path(bif_dir, paste0(c, '.rds')))
        rm(npa_list)
        rm(bif)
    })
})
