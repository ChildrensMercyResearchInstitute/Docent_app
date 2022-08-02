FROM rocker/r-base:latest
LABEL maintainer="manish <mkumar1@cmh.edu>"
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
    shiny \
    forecast \
    jsonlite \
    ggplot2 \
    htmltools \
    plotly \  
    shinyjs \
    shinyWidgets \
    shinycssloaders \
    magrittr \
    dplyr \
    tidyr \
    tibble \
    glue \
    parallel 

RUN install2.r --error --skipinstalled \
    yaml \
    RSQLite \
    readr \ 
    optparse \
    tibble 
 
RUN R -e "install.packages('VariantAnnotation',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"


RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site
RUN addgroup --system app \
    && adduser --system --ingroup app app
WORKDIR /home/app
COPY docent .
COPY log .
RUN chown app:app -R /home/app
USER app
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/home/app')"]
