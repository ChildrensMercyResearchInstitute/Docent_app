FROM rocker/r-base:latest
LABEL maintainer="USER <user@example.com>"
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*
#FROM rocker/shiny:4.0.5
# Install system requirements for index.R as needed
#RUN apt-get update && apt-get install -y \
#    --no-install-recommends \
#    git-core \
#    libssl-dev \
#    libcurl4-gnutls-dev \
#    curl \
#    libsodium-dev \
#    libxml2-dev \
#    libicu-dev \
#    libsqlite3-dev \ 
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/*
#ENV _R_SHLIB_STRIP_=true
#COPY Rprofile.site /etc/R
RUN install2.r --error --skipinstalled \
    shiny \
    forecast \
    jsonlite \
    ggplot2 \
    htmltools \
    plotly  

RUN R -e "install.packages('shiny',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyjs',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyWidgets',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinycssloaders',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magrittr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tibble',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('mongolite',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('glue',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('parallel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('yaml',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('RSQLite',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('VariantAnnotation',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magrittr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('raedr',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tibble',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('mongolite',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site
RUN addgroup --system app \
    && adduser --system --ingroup app app
WORKDIR /home/app
COPY docent_poster .
COPY log .
RUN chown app:app -R /home/app
USER app
EXPOSE 3838
#CMD ["sh"]

CMD ["R", "-e", "shiny::runApp('/home/app')"]

#COPY docent_poster /srv/shiny-server/
#COPY log /srv/shiny-server/
#USER shiny
#EXPOSE 3838
#CMD ["/usr/bin/shiny-server"]
