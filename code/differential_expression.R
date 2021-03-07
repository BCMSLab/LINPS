# perform differential expression

# load required libraries
library(SummarizedExperiment)
library(limma)
library(tidyverse)

# load data
perturbs <- c("trt_sh", "trt_oe",  "trt_oe.mut", "trt_cp")
controls <- c(rep('EMPTY_VECTOR', 3), 'DMSO')

res_dir <- 'results/'

map2(perturbs, controls,
     function(trt, ctr) {
         trt_dir <- paste0(res_dir, trt)
         dir.create(trt_dir)
         
         res <- paste0(trt_dir, '/diff_expr')
         dir.create(res)
         message(paste('Directory', res, 'created.'))
         
         d <- paste0(res_dir, trt, '/exprs/')
         cell <- str_split(list.files(d),
                           '\\.',
                           simplify = TRUE)[, 1]
         
         map(cell, function(c) {
             se_ctr <- read_rds(paste0(res_dir, ctr, '/exprs/', c, '.rds'))
             se_trt <- read_rds(paste0(res_dir, trt, '/exprs/', c, '.rds'))
             
             se <- cbind(se_trt, se_ctr)
             
             message(paste('SE object', trt, 'in', c, 'is loaded.'))
             rm(se_ctr)
             rm(se_trt)
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
             rm(se)
             
             message(paste('Done. Differential expression of', trt, 'in', c))
         })
     })
