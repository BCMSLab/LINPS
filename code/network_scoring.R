# loading required libraries
library(tidyverse)
library(reshape2)
library(NPA)
library(NPAModels)
library(BiocParallel)

# load models
models <- load_models('Hs')

# create a results directory
res_dir <- file.path('results')
diff_expr <- file.path(res_dir, 'diff_expr')
if (!file.exists(diff_expr)) dir.create(diff_expr)

npa_list <- file.path(res_dir, 'npa_list')
if (!file.exists(npa_list)) dir.create(npa_list)
bif <- file.path(res_dir, 'bif')
if (!file.exists(bif)) dir.create(bif)

perturbs <- list.files(diff_expr)

map(perturbs, function(x) {
    message(paste('Scoring', x, 'perturbations.'))
    trt_dir <- file.path(diff_expr, x)
    
    # create directories
    npa_list_dir <- file.path(npa_list, x)
    if (!dir.exists(npa_list_dir)) dir.create(npa_list_dir)
    
    bif_dir <- file.path(res_dir, 'bif', x)
    if (!dir.exists(bif_dir)) dir.create(bif_dir)
    
    if (length(dir(trt_dir)) != 0) {
        cell_lines <- str_split(list.files(trt_dir),
                                '\\.',
                                simplify = TRUE)[, 1]
    
        map(cell_lines, function(cell) {
            # read in differential expression table
            diff_expr <- read_tsv(file.path(trt_dir, paste0(cell, '.tsv')))
            
            # make groups (for NPA)
            groups <- diff_expr %>%
                dplyr::select(nodeLabel = ID,
                              t,
                              foldChange = logFC) %>%
                as.data.frame() %>%
                with(split(., diff_expr$pert_iname)) %>%
                map(function(x) {
                    filter(x, !duplicated(x$nodeLabel))
                })
            
            # check names are correct
            all(map(groups, ~all(sum(duplicated(toupper(.x))))) == 0)
            
            # clean
            rm(diff_expr)
            
            # compute npa
            npa_list <- compute_npa_list(
                groups,
                models,
                verbose = TRUE)
            
            # write to disk
            write_rds(npa_list, file.path(npa_list_dir, paste0(cell, '.rds')))
            
            # clean
            rm(groups)
            
            # compute bif
            bif <- get_bif(npa_list)
            
            # write to disk
            write_rds(bif, file.path(bif_dir, paste0(cell, '.rds')))
            
            # clean
            rm(npa_list)
            rm(bif)
        })    
    }
})
