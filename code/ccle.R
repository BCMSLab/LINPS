# getGEOSuppFiles('GSE36133', fetch_files = FALSE)
# aria2c -x 8 -s 8 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36133/suppl//GSE36133_RAW.tar

ccle <- getGEO('GSE36133', destdir = 'data/')[[1]] 

pData(ccle) %>%
    as_tibble() %>%
    select(cell_id = title, 
           tissue = `primary site:ch1`) %>%
    unique() %>%
    mutate(cell_id = toupper(str_replace(cell_id, '\\-| ', '')),
           tissue = str_to_title(tissue)) %>%
    write_tsv('results/cell_tissue.tsv')
