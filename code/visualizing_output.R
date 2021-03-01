library(tidyverse)

con <- DBI::dbConnect(RSQLite::SQLite(),
                      dbname = "results/LINPS.sqlite")

dbListTables(con)

# the relative biological impact factor
tbl(con, 'bif') %>%
    filter(cell_id == 'MCF7',
           pert_type == 'trt_cp',
           type == 'rbif') %>%
    collect() %>%
    mutate(pert_iname = fct_reorder(pert_iname, desc(BIF))) %>%
    top_n(BIF, n = Inf) %>%
    ggplot(aes(x = pert_iname, y = BIF)) +
    geom_col()

# the biological impact factor per network
families <- c('CPR', 'IPN')

tbl(con, 'bif') %>%
    filter(cell_id == 'MCF7',
           pert_type == 'trt_cp',
           type == 'r2') %>%
    collect() %>%
    gather(family, r2, families) %>%
    mutate(pert_iname = fct_reorder(pert_iname, desc(r2))) %>%
    top_n(r2, n = Inf) %>%
    ggplot(aes(x = pert_iname, y = r2, fill = family)) +
    geom_col()

# npa for a given network
tbl(con, 'networks') %>%
    filter(cell_id == 'MCF7',
           pert_type == 'trt_cp',
           network == 'Cell_Cycle') %>%
    collect() %>%
    mutate(pert_iname = fct_reorder(pert_iname, desc(coefficients))) %>%
    top_n(coefficients, n = Inf) %>%
    ggplot(aes(x = pert_iname, y = coefficients)) +
    geom_col() +
    geom_errorbar(aes(ymin = ci.down, ymax = ci.up))

# nodes
tbl(con, 'nodes') %>%
    filter(cell_id == 'MCF7',
           pert_type == 'trt_cp',
           network == 'Cell_Cycle',
           nodes.coefficients.pvalue < .01) %>%
    collect() %>%
    mutate(node = fct_reorder(node, desc(abs(nodes.coefficients)))) %>%
    group_by(pert_iname) %>%
    top_n(nodes.coefficients, n = 5) %>%
    ggplot(aes(x = node,
               y = nodes.coefficients,
               ymin = nodes.coefficients.ci.down,
               ymax = nodes.coefficients.ci.up)) +
    geom_col() +
    geom_errorbar() +
    facet_wrap(~pert_iname, scales = 'free_x')

dbDisconnect(con)

# networks
# load lincs info
db_fl <- "LINPS/results/LINPS.sqlite"
# db_fl <- "www/LINPS.sqlite"
file.exists(db_fl)

con <- dbConnect(SQLite(), dbname = db_fl)

cell_lines <- tbl(con, 'perturbations') %>%
    select(cell_id) %>%
    collect() %>%
    unique() %>%
    mutate(tissue = c('Skin', 'Breast')) # TODO: to be removed

perturbations <- tbl(con, 'perturbations') %>%
    collect() %>%
    unique()

models <- tbl(con, 'models') %>%
    collect() %>%
    filter(model %in% c('Cell_Cycle',
                        'Epithelial_Innate_Immune_Activation')) # TODO: to be removed

graphs <- map(models$graph, unserialize)
names(graphs) <- models$model

dbDisconnect(con)

g <- graphs[[2]]

n <- as_data_frame(g, 'vertices') %>%
    select(id = name)

e <- as_data_frame(g, 'edges') %>%
    as_tibble() %>%
    mutate(color = c('green', 'red')[factor(e$weight)])

visNetwork(nodes = n,
           edges = e) %>%
    visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 2))) %>%
    visOptions(nodesIdSelection = TRUE,
               highlightNearest = TRUE) %>%
    visLayout(randomSeed = 123)
