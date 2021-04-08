change_ids <- function(se, rownames, org) {
    id_symbol <- AnnotationDbi::select(org,
                                       rownames(se),
                                       rownames,
                                       'ENTREZID')
    
    # check rownames in correspond to entrez ids
    stopifnot(all(rownames(se) == id_symbol$ENTREZID))
    
    # use symbols as rownames
    rownames(se) <- id_symbol$SYMBOL
    
    # return 
    return(se)
}

subset_sl <- function(sl, col.ix, n, rownames = 'SYMBOL', org, file) {
    # check col.ix is not 0
    stopifnot(length(col.ix) > 0)
    stopifnot(rownames %in% AnnotationDbi::columns(org))
    
    # sample col.ix when required
    if (!missing(n)) {
        col.ix <- sample(col.ix, min(length(col.ix), n))
    }
    
    # subset sl object
    se <- as(sl[, col.ix], "SummarizedExperiment")
    
    # map ids to symbols
    if (!missing(org)) {
        se <- change_ids(se, 'SYMBOL', org)
    }
    
    # remove rows with na in the row names
    se <- se[!is.na(rownames(se)),]
    
    # return
    if (!missing(file)) {
        readr::write_rds(se, file)
    } else {
        return(se)
    }
}

differential_expression <- function(se, control_name, file) {
    # make factor
    se$pert_iname <- relevel(factor(se$pert_iname),
                             ref = control_name)
    
    # remove missing values
    se <- se[!is.na(rownames(se)),]
    
    # TODO: a more sophisticated model needs to be used
    mod <- model.matrix(~pert_iname, colData(se))
    colnames(mod) <- str_replace(colnames(mod), 'pert_iname', '')
    
    # apply lmFit
    fit <- lmFit(assay(se), mod)
    fit <- eBayes(fit)
    
    # get toptable    
    df <- map_df(colnames(mod)[-1],
                 ~topTable(fit,
                           number = Inf,
                           genelist = rownames(se),
                           coef = .x) %>%
                     mutate(pert_iname = .x)) %>%
        mutate(cell_id = unique(se$cell_id))
    
    # return
    if (!missing(file)) {
        write_tsv(df, file)
    } else {
        return(df)
    }
}