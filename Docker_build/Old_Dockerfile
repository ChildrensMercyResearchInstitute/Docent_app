FROM ubuntu:18.04
ENV R_BASE_VERSION 4.1.0
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
		r-base

MAINTAINER Manish Kumar
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

RUN R -e "install.packages('VariantAnnotation',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magrittr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('raedr',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tibble',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN apt-get install -y sqlite3 libsqlite3-dev


CMD ["/bin/sh"]
