########### Dockerfile, for use with scilifelab serve

FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Command to install standard R packages from CRAN; enter the list of required packages for your app here
RUN Rscript -e 'install.packages(c("shiny","tidyverse","plotly","Cairo","shinyjs","patchwork","sqldf","reshape2","cowplot", "umap","stringr"))'


# Command to install packages from Bioconductor; enter the list of required Bioconductor packages for your app here
#RUN Rscript -e 'BiocManager::install(c("Biostrings"),ask = F)'

RUN rm -rf /srv/shiny-server/*
#COPY /app/ /srv/shiny-server/
COPY / /srv/shiny-server/

USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
