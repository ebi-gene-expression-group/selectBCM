FROM r-base


RUN mkdir /home/analysis


RUN apt-get -y update && apt-get install -y \
    libssl-dev \
    libcairo2-dev \
    libjpeg-dev \
    cmake \
    libgif-dev \
    build-essential \
    libxml2-dev\
    libfontconfig1-dev\
    curl\
    libcurl4-openssl-dev\
    libpng-dev\
    default-jdk \
    r-cran-rjava \
    git



   RUN R -e "system('apt-get install -y libxml2-dev')"
   RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
   RUN R -e 'options(download.file.method = "libcurl")'
   RUN R -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")'


 ENV RENV_VERSION 0.16.0

 RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

  ## add installation for scBatch
  RUN R -e 'remotes::install_github("tengfei-emory/scBatch", ref="master")'

   ## Install SelectBCM from github
  RUN R -e 'remotes::install_github("ebi-gene-expression-group/selectBCM", ref="master")'


 RUN  apt-get -y autoclean && \
   rm -rf /var/lib/apt/lists/* && \
   rm -rf /tmp/* \
   && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
