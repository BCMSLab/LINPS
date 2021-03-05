# download the gene expression signatures dataset
# make a metadata table
# make a test object ++ to be removed

# download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_inst_info.txt.gz',
#               destfile = 'data/GSE92742_Broad_LINCS_inst_info.txt.gz')
# aria2c -x 8 -s 8 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz

# load libraries
library(tidyverse)
library(org.Hs.eg.db)
library(slinky)

# set an sl instance
user_key <- httr::content(httr::GET(
    "https://api.clue.io/temp_api_key"),
    as = "parsed")$user_key

data_directory <- '~/workingon/anti_metastatic3/data' ## to be changed

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

# TODO: select a test set ++ to be removed
# extract metadata table
md <- metadata(sl)

cell_lines <- head(pull(md, cell_id) %>% unique(), 10)
control_n <- 100
perturbations_n <- 100

# subset the data

perturbs <- c("trt_cp", "trt_sh", "trt_oe",  "trt_oe.mut")
controls <- c('DMSO', rep('EMPTY_VECTOR', 3))

map2(perturbs, controls,
     function(trt, ctr) {
         # subset perturbations
         trt_ind <- head(filter(md, pert_type == trt) %>% pull(pert_iname) %>% unique(), perturbations_n)
         
         # write md file
         md_file <- paste0('data/expression_metadata_', trt, '.tsv')
         md %>%
             dplyr::select(starts_with('pert'), cell_id, -pert_id) %>%
             unique() %>%
             filter(pert_iname %in% trt_ind,
                    cell_id %in% cell_lines) %>%
             write_tsv(md_file)
         
         # subset sl object
         col.ix <- which(md$cell_id %in% cell_lines & md$pert_iname %in% trt_ind)
         set.seed(123)
         col.ix2 <- sample(which(md$cell_id %in% cell_lines & md$pert_iname == ctr))

         # se <- as(sl[, c(col.ix, col.ix2[1:control_n]) ], "SummarizedExperiment")
         se <- as(sl[, c(col.ix, col.ix2) ], "SummarizedExperiment")
         
         # map ids to symbols
         id_symbol <- select(org.Hs.eg.db,
                             rownames(se),
                             'SYMBOL',
                             'ENTREZID')
         all(rownames(se) == id_symbol$ENTREZID)
         rownames(se) <- id_symbol$SYMBOL

         se <- se[!is.na(rownames(se)),]

         # write to file
         write_rds(se, paste0('data/', trt, '.rds'))
     })

# merge the metadatafiles
list.files('data',
           pattern = 'expression_metadata_trt*', full.names = TRUE) %>%
    map_df(function(x) {
        df <- read_tsv(x) %>% mutate(pert_dose_unit = as.character(pert_dose_unit))
        unlink(x)
        df
    }) %>%
    write_tsv('data/expression_metadata.tsv')
