# load the required libraries
library(tidyverse)
library(reshape2)
library(DBI)

# make the database file
if(file.exists('results/LINPS.sqlite')) unlink('results/LINPS.sqlite')
con <- dbConnect(RSQLite::SQLite(), dbname = "results/LINPS.sqlite")

# differential expression
# list.dirs('results', recursive = FALSE, full.names = FALSE) %>%
#     map(function(trt) {
#         d1 <- paste('results', trt, sep = '/')
#         fls <- list.files(d1,
#                           pattern = '*.tsv', 
#                           recursive = TRUE,
#                           full.names = TRUE)
#         fls %>%
#             map(function(f) {
#                 df <- read_tsv(f)
#                 
#                 # TODO: should be added in diff_expr
#                 df <- mutate(df, pert_type = trt)
#                 
#                 dbWriteTable(con,
#                              value = df,
#                              name = 'diff_expr',
#                              overwrite = FALSE,
#                              append = TRUE)
#             })
#     })

# biological impact factor
list.dirs('results', recursive = FALSE, full.names = FALSE) %>%
    map(function(trt) {
        d1 <- paste('results', trt, 'bif', sep = '/')
        fls <- list.files(d1,
                          pattern = '*.rds', 
                          recursive = TRUE,
                          full.names = TRUE)
        fls %>%
            map(function(f) {
                fl <- str_split(f,
                                '/|\\.',
                                simplify = TRUE)
                bif <- read_rds(f)
                
                # TODO: extract other tables from BIF
                # TODO: extract other items in the object
                df <- map_df(bif$get_data()$BIF[c('rbif', 'r2')],
                             function(x) melt(x) %>% setNames(c('pert_iname', 'stat', 'value')),
                             .id = 'type') %>%
                    mutate(pert_type = fl[, 2],
                           cell_id = fl[, 4]) %>%
                    dplyr::select(pert_type, pert_iname, cell_id,
                                  everything()) %>%
                    spread(stat, value)
                dbWriteTable(con,
                             value = df,
                             name = 'bif',
                             overwrite = FALSE,
                             append = TRUE)
            })
    })

list.dirs('results', recursive = FALSE, full.names = FALSE) %>%
    map(function(trt) {
        d1 <- paste('results', trt, 'npa_lists', sep = '/')
        fls <- list.files(d1,
                          pattern = '*.rds', 
                          recursive = TRUE,
                          full.names = TRUE)
        fls
        fls %>%
            map(function(f) {
                fl <- str_split(f,
                                '/|\\.',
                                simplify = TRUE)
                npa_list <- read_rds(f)
                
                # TODO: extract other items in the object
                # network info
                net <- c('coefficients', 'coefficients.var', 'ci.up', 'ci.down')
                
                networks <- npa_list$get_data() %>%
                    map_df(function(n) {
                        bind_rows(n[net]) %>%
                            mutate(pert_iname = names(n$coefficients))
                    }, .id = 'network') %>%
                    separate(network,
                             into = c('family', 'network'),
                             sep = ' / ') %>%
                    mutate(pert_type = fl[, 2],
                           cell_id = fl[, 4]) %>%
                    dplyr::select(pert_type, pert_iname, cell_id,
                                  family, network,
                                  everything())
                dbWriteTable(con,
                             value = networks,
                             name = 'networks',
                             overwrite = FALSE,
                             append = TRUE)
                
                # nodes info
                nds <- c("nodes.coefficients", "nodes.coefficients.ci.up",
                         "nodes.coefficients.ci.down", "nodes.coefficients.pvalue" )
                nodes <- npa_list$get_data() %>%
                    map_df(function(n) {
                        map_df(n[nds], melt, .id = 'stat') %>% as_tibble()
                    }, .id = 'network') %>%
                    separate(network,
                             into = c('family', 'network'),
                             sep = ' / ') %>%
                    mutate(pert_type = fl[, 2],
                           cell_id = fl[, 4]) %>%
                    dplyr::select(pert_type, pert_iname = Var2, cell_id,
                                  family, network, node = Var1,
                                  everything()) %>%
                    spread(stat, value)
                dbWriteTable(con,
                             value = nodes,
                             name = 'nodes',
                             overwrite = FALSE,
                             append = TRUE)
            })
    })

# the network models
models <- read_tsv('data/models_metadata.tsv') %>%
    mutate(version = str_replace_all(version, '\\.', '__'),
           graph = paste(species, family, model, version, sep = '__'))

models$graph <- models$graph %>%
    map(function(x) {
        fl <- file.path('data/NPAModels/data', paste0(x, '.rda'))
        
        env <- new.env()
        load(fl, envir = env)
        mod <- env[[ls(env)[[1]]]]
        serialize(mod$g, NULL)
    })

dbWriteTable(con,
             value = models,
             name = 'models',
             overwrite = FALSE,
             append = TRUE)

# the perturbation metadata
perturbations <- read_tsv('data/expression_metadata.tsv')
dbWriteTable(con,
             value = perturbations,
             name = 'perturbations',
             overwrite = FALSE,
             append = TRUE)

# disconnect
dbDisconnect(con)
