# Dockerfile for the SelectBCM v2.0.0
# Written and maintain by Madhulika Mishra(madhulika@ebi.ac.uk) 
FROM r-base


RUN mkdir /home/analysis


RUN apt-get -y update && apt-get install -y \
    libssl-dev \
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
    libxt-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    libx11-dev \
    r-cran-nloptr \
    libcairo2-dev \
    git



   RUN R -e "system('apt-get install -y libxml2-dev')"
   RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
   RUN R -e 'options(download.file.method = "libcurl")'
   RUN R -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")'


  ## add installation for scBatch
  RUN R --no-echo --no-restore --no-save -e 'remotes::install_github("tengfei-emory/scBatch", ref="master")'

   ## Install SelectBCM from github
  RUN R  --no-echo --no-restore --no-save  -e 'remotes::install_github("ebi-gene-expression-group/selectBCM", ref="master",build = FALSE)'


 RUN  apt-get -y autoclean && \
   rm -rf /var/lib/apt/lists/* && \
   rm -rf /tmp/* \
   && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

CMD [ "R" ]
