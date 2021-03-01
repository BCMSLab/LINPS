# get and preprocess the causal networks (NPAModels)
# make a metadata table

# TODO: more causal networks can be added

# load required libraries
library(tidyverse)
library(NPAModels)

# install the NPAModles
devtools::install('data/NPAModels/')

# write a metadata table
tibble(file = list.files('data/NPAModels/data/', 'Hs')) %>%
    separate(file,
             c('species', 'family', 'model',
               paste0('version_', 1:3)), 
             '__|\\.') %>%
    unite('version', starts_with('version'), sep = '.') %>%
    write_tsv('data/models_metadata.tsv')

# preprocess the networks
preprocessNetworks()

