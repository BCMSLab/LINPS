# perform differential expression

# load required libraries
library(SummarizedExperiment)
library(limma)
library(tidyverse)

# load data
perturbs <- c("trt_cp", "trt_sh", "trt_oe",  "trt_oe.mut")
controls <- c('DMSO', rep('EMPTY_VECTOR', 3))

map2(perturbs, controls,
     function(trt, ctr) {
         trt_dir <- paste0('results/', trt)
         dir.create(trt_dir)
         
         res <- paste0(trt_dir, '/diff_expr')
         dir.create(res)
         message(paste('Directory', res, 'created.'))
         
         trt_se <- read_rds(file.path('data', paste0(trt, '.rds')))
         message(paste('SE object', trt, 'is loaded.'))
         
         map(unique(trt_se$cell_id), function(c) {
             se <- trt_se[, trt_se$cell_id == c]
             se$pert_iname <- relevel(factor(se$pert_iname), ref = ctr)
             se <- se[!is.na(rownames(se)),]
             
             # TODO: a more sophisticated model needs to be used
             mod <- model.matrix(~pert_iname, colData(se))
             colnames(mod) <- str_replace(colnames(mod), 'pert_iname', '')
             
             fit <- lmFit(assay(se), mod)
             fit <- eBayes(fit)
             
             map_df(colnames(mod)[-1],
                 ~topTable(fit,
                           number = Inf,
                           genelist = rownames(se),
                           coef = .x) %>%
                     mutate(pert_iname = .x)) %>%
                 mutate(cell_id = c) %>%
                 write_tsv(file.path(res, paste0(c, '.tsv')))
             
             rm(fit)
             
             message(paste('Done. Differential expression of', trt, 'in', c))
         })
         rm(trt_se)
     })
