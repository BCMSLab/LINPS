# load the required libraries
library(tidyverse)
library(reshape2)
library(DBI)

# make the database file
res_dir <- 'results/'
fl <- file.path(res_dir, 'LINPS.sqlite')

if(file.exists(fl)) unlink(fl)

con <- dbConnect(RSQLite::SQLite(), dbname = fl)

# TODO: add significance (pvalues)
# biological impact factor
bif_dir <- file.path('results', 'bif')
npa_list_dir <- file.path('results', 'npa_list')
sim_dir <- file.path(res_dir, 'similarity')

trt <- list.files(bif_dir)

trt %>%
    map(function(t) {
        
        fls <- list.files(file.path(bif_dir, t),
                          pattern = '*.rds', 
                          full.names = TRUE)
        
        map(fls, function(f) {
            # extract info from file name
            fl <- str_split(f,
                            '/|\\.rds',
                            simplify = TRUE)
            
            # load object
            bif <- read_rds(f)
            
            # extract and write to the db
            tibble(pert_type = str_replace_all(fl[, 3], '\\.', ' '),
                   pert_iname = rownames(NPA::as.matrix(bif)),
                   cell_id = fl[, 4],
                   RBIF = NPA::as.matrix(bif, type = 'coefficients')[, 'RBIF']) %>%
                filter(RBIF > 0) %>%
                dbWriteTable(con,
                             value = .,
                             name = 'rbif',
                             overwrite = FALSE,
                             append = TRUE)
            
            # extract and write to the db
            NPA::as.matrix(bif) %>%
                as.data.frame() %>%
                dplyr::select(-BIF) %>%
                filter(rowSums(.) > 0) %>%
                rownames_to_column('pert_iname') %>%
                as_tibble() %>%
                mutate(pert_type = str_replace_all(fl[, 3], '\\.', ' '),
                       cell_id = fl[, 4]) %>%
                dplyr::select(pert_type, pert_iname, cell_id, everything()) %>%
                dbWriteTable(con,
                             value = .,
                             name = 'bif',
                             overwrite = FALSE,
                             append = TRUE)
        })
    })

trt %>%
    map(function(t) {
        fls <- list.files(file.path(npa_list_dir, t),
                          pattern = '*.rds', 
                          full.names = TRUE)
        
        fls[lengths(fls) > 0] %>%
            map(function(f) {
                fl <- str_split(f,
                                '/|\\.rds',
                                simplify = TRUE)
                
                npa_list <- read_rds(f)
                
                # network info
                nets <- NPA::networks(npa_list)
                
                map(seq_along(nets), function(x) {
                    npa_sub <- NPA::subset(npa_list, x)
                    tibble(
                        pert_iname = NPA::comparisons(npa_sub),
                        coefficients = NPA::coefficients(npa_sub),
                        down = NPA::conf.int(npa_sub)[, 'down'],
                        up = NPA::conf.int(npa_sub)[, 'up']
                    ) %>%
                        mutate(pert_type = str_replace_all(fl[, 3], '\\.', ' '),
                               cell_id = fl[, 4],
                               network = nets[x]) %>%
                        separate(network,
                                 into = c('family', 'network'),
                                 sep = ' / ') %>%
                        dplyr::select(pert_type, pert_iname, cell_id, 
                                      family, network, everything()) %>%
                        dbWriteTable(con,
                                     value = .,
                                     name = 'networks',
                                     overwrite = FALSE,
                                     append = TRUE)
                    
                    # nodes
                    # extract coefficients
                    coefficients <-  NPA::coefficients(npa_sub, type = 'nodes') %>%
                        melt() %>%
                        setNames(c('node', 'pert_iname', 'coefficient')) %>%
                        mutate_at(vars('node', 'pert_iname'), as.character)
                    
                    # extract conf.int
                    conf.int <- NPA::conf.int(npa_sub, type = 'nodes') %>%
                        melt() %>%
                        separate(Var2, into = c('Var2', 'dir'), sep = ' \\(') %>%
                        mutate(dir = str_remove_all(dir, '\\)')) %>%
                        spread(dir, value) %>%
                        setNames(c('node', 'pert_iname', 'down', 'up')) %>%
                        mutate_at(vars('node', 'pert_iname'), as.character)
                    
                    # join, format and wirte to db
                    inner_join(coefficients, conf.int) %>%
                        mutate(pert_type = fl[, 3],
                               pert_type = str_replace_all(pert_type, '\\.', ' '),
                               cell_id = fl[, 4],
                               network = nets[x]) %>%
                        separate(network,
                                 into = c('family', 'network'),
                                 sep = ' / ') %>%
                        dplyr::select(pert_type, pert_iname, cell_id, 
                                      family, network, node, everything()) %>%
                        dbWriteTable(con,
                                     value = .,
                                     name = 'nodes',
                                     overwrite = FALSE,
                                     append = TRUE)
                }) 
            })
    })

# similarity
list.files('results/similarity', recursive = TRUE, full.names = TRUE) %>%
    map(function(f) {
        df <- read_tsv(f) %>% 
            filter(PVAL < .05)
        dbWriteTable(con,
                     value = df,
                     name = 'similarity',
                     overwrite = FALSE,
                     append = TRUE)
    })

# the network models
models <- read_tsv(file.path(res_dir, 'models_metadata.tsv')) %>%
    mutate(version = str_replace_all(version, '\\.', '__'),
           graph = paste(species, family, model, version, sep = '__'))

models$graph <- models$graph %>%
    map(function(x) {
        fl <- file.path('tools/NPAModels/data', paste0(x, '.rda'))
        
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
perturbations <- read_tsv(file.path(res_dir, 'expression_metadata.tsv'))

dbWriteTable(con,
             value = perturbations,
             name = 'perturbations',
             overwrite = TRUE)

# disconnect
dbDisconnect(con)
