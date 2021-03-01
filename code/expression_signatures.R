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
trt <- head(filter(md, pert_type == 'trt_cp') %>% pull(pert_iname) %>% unique(), 10)
sh <- head(filter(md, pert_type == 'trt_sh') %>% pull(pert_iname) %>% unique(), 10)
cell_lines <- head(pull(md, cell_id) %>% unique(), 5)

# trt <- c('rutin', 'albendazole', 'tamoxifen', 'alrestatin', 'altizide', 'altrenogest')
# sh <- c('ADPGK', 'AK3', 'AKR1C2', 'BAD', 'CDCA7L')
# cell_lines <- c('MCF7', 'A375')

md %>%
    dplyr::select(starts_with('pert'), cell_id, -pert_id) %>%
    unique() %>%
    filter(pert_iname %in% c(trt, sh),
           cell_id %in% cell_lines) %>%
    write_tsv('data/expression_metadata.tsv')

col.ix <- which(metadata(sl)$cell_id %in% cell_lines & metadata(sl)$pert_iname %in% trt)
col.ix2 <- which(metadata(sl)$cell_id %in% cell_lines & metadata(sl)$pert_iname == 'UnTrt')

# trt cp
trt_cp <- as(sl[, c(col.ix, col.ix2[1:100]) ], "SummarizedExperiment")

id_symbol <- select(org.Hs.eg.db,
                    rownames(trt_cp),
                    'SYMBOL',
                    'ENTREZID')
all(rownames(trt_cp) == id_symbol$ENTREZID)
rownames(trt_cp) <- id_symbol$SYMBOL

trt_cp <- trt_cp[!is.na(rownames(trt_cp)),]

write_rds(trt_cp, 'data/trt_cp.rds')

col.ix <- which(metadata(sl)$cell_id %in% cell_lines & metadata(sl)$pert_iname %in% sh)
col.ix2 <- which(metadata(sl)$cell_id %in% cell_lines & metadata(sl)$pert_iname == 'EMPTY_VECTOR')

# trt sh
trt_sh <- as(sl[, c(col.ix, col.ix2[1:100]) ], "SummarizedExperiment")

id_symbol <- select(org.Hs.eg.db,
                    rownames(trt_sh),
                    'SYMBOL',
                    'ENTREZID')
all(rownames(trt_sh) == id_symbol$ENTREZID)
rownames(trt_sh) <- id_symbol$SYMBOL

write_rds(trt_sh, 'data/trt_sh.rds')
