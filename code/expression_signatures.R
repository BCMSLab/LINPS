# load libraries
library(tidyverse)
library(org.Hs.eg.db)
library(slinky)
source('R/functions.R')

# set an sl instance
user_key <- httr::content(httr::GET(
    "https://api.clue.io/temp_api_key"),
    as = "parsed")$user_key

# TODO: to be changed
data_directory <- 'data' 

file_gctx <- file.path(
    data_directory,
    'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx')
file.exists(file_gctx)

file_info <- file.path(
    data_directory,
    'GSE92742_Broad_LINCS_inst_info.txt.gz')
file.exists(file_info)

sl <- Slinky(
    user_key,
    file_gctx,
    file_info)

# extract metadata table
compinfo <- read_tsv('data/compoundinfo_beta.txt')
tissue_cell <- read_tsv('data/ccle_tisse_cells.tsv')
    
md <- metadata(sl)
md_sub <- md %>% 
    mutate(rowid = row_number()) %>%
    as_tibble() %>%
    left_join(compinfo) %>%
    left_join(tissue_cell) %>%
    mutate(tissue = ifelse(is.na(tissue), 'other', tissue)) %>%
    filter(pert_type == 'trt_cp', !is.na(moa)) %>%
    group_by(moa, cell_id, pert_type, pert_iname) %>%
    summarise(idx = list(unique(rowid)),
              n = length(unique(unlist(idx)))) %>%
    ungroup() %>%
    filter(n > 3) %>%
    unique()

# create results directories
res_dir <- file.path('results')
if (!dir.exists(res_dir)) dir.create(res_dir)

md[unique(unlist(md_sub$idx)),] %>%
    dplyr::select(-pert_type) %>%
    left_join(tissue_cell) %>%
    left_join(dplyr::select(compinfo, pert_id, pert_type = moa)) %>%
    as_tibble() %>%
    write_tsv(file.path(res_dir, 'expression_metadata.tsv'))

# subset and write expression sets
exprs_dir <- file.path(res_dir, 'exprs')

if (!dir.exists(exprs_dir)) dir.create(exprs_dir)

md_sub %>%
    group_split(moa) %>%
    map(function(x) {
        d <- file.path(exprs_dir, make.names(unique(x$moa)))
        if (!file.exists(d)) dir.create(d)
        
        cells <- unique(x$cell_id)
        
        map(cells, function(c) {
            fl <- file.path(d, paste0(c, '.rds'))
            ids <- filter(x, cell_id == c) %>% pull(idx) %>% unlist()
            subset_sl(
                sl,
                ids,
                rownames = 'SYMBOL',
                org = org.Hs.eg.db,
                file = fl
            )

            # message
            message(paste('Done.', unique(x$moa), 'in', c))
        })
    })

# controls
controls_sub <- md %>% 
    mutate(rowid = row_number()) %>%
    as_tibble() %>%
    filter(pert_type == 'ctl_vehicle', pert_iname == 'DMSO') %>%
    group_by(cell_id) %>%
    summarise(idx = list(head(unique(rowid), 100)))

controls_dir <- file.path(exprs_dir, 'DMSO')
if (!dir.exists(controls_dir)) dir.create(controls_dir)

controls_sub %>%
    group_split(cell_id) %>%
    map(function(c) {
        fl <- file.path(controls_dir, paste0(c$cell_id, '.rds'))
        ids <-  unlist(c$idx)
        
        subset_sl(
            sl,
            ids,
            rownames = 'SYMBOL',
            org = org.Hs.eg.db,
            file = fl
        )
        
        # message
        message(paste('Done. DMSO', 'in', c$cell_id))
    })
