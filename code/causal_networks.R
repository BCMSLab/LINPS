# get and preprocess the causal networks (NPAModels)
# make a metadata table

# load required libraries
library(tidyverse)
library(NPAModels)

# write a metadata table
tibble(file = list.files('tools/NPAModels/data/', 'Hs')) %>%
    separate(file,
             c('species', 'family', 'model',
               paste0('version_', 1:3)), 
             '__|\\.') %>%
    unite('version', starts_with('version'), sep = '.') %>%
    write_tsv('data/models_metadata.tsv')

# preprocess the networks
preprocessNetworks()
