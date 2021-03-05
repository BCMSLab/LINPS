# get shiny serves plus tidyverse packages image
FROM bioconductor/bioconductor_docker:latest

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev 
  
# install R packages required 
## required by bioc

RUN R -e "BiocManager::install(c('SummarizedExperiment', 'slinky', 'org.Hs.eg.db', 'limma'))"

## other tidy pkgs
RUN R -e "install.packages(c('tidyverse', 'reshape2'))"

## required by the app
RUN R -e "install.packages(c('shiny', 'shinyBS', 'DT'))"
RUN R -e "install.packages(c('RSQLite', 'DBI'))"
RUN R -e "install.packages(c('plotly', 'igraph', 'visNetwork', 'Hmisc'))"

## required from github
RUN R -e "remotes::install_github('MahShaaban/NPAModels@hsonly')"
RUN R -e "remotes::install_github('philipmorrisintl/NPA')"
