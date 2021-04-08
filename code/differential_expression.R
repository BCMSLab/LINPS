# perform differential expression

# load required libraries
library(SummarizedExperiment)
library(limma)
library(tidyverse)
source('R/functions.R')

# create directory
res_dir <- file.path('results')
exprs_dir <- file.path(res_dir, 'exprs')
diff_expr <- file.path(res_dir, 'diff_expr')
if (!file.exists(diff_expr)) dir.create(diff_expr)

trt <- list.files(exprs_dir)
trt <- trt[!grepl('DMSO', trt)]

map(trt, function(t) {
    cell_lines <- str_split(list.files(file.path(exprs_dir, t)),
                            '\\.',
                            simplify = TRUE)[, 1]
    trt_deg <- file.path(diff_expr, t)
    if (!file.exists(trt_deg)) dir.create(trt_deg)
    
    map(cell_lines, function(c) {
        ctr_file <- file.path(exprs_dir, 'DMSO', paste0(c, '.rds'))
        trt_file <- file.path(exprs_dir, t, paste0(c, '.rds'))
        
        # read expression
        se_ctr <- read_rds(ctr_file)
        se_trt <- read_rds(trt_file)
        
        # combine
        se <- cbind(se_trt, se_ctr)
        
        # clean
        rm(se_ctr)
        rm(se_trt)
        
        # write top tables
        fl <- file.path(trt_deg, paste0(c, '.tsv'))
        
        # get deg
        differential_expression(se, 'DMSO', fl)
        
        # # clean
        # rm(se)
        # 
        message(paste('Done. Differential expression of', t, 'in', c))
    })
})
