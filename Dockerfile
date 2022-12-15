FROM rocker/r-ver:4.2.2
RUN apt-get update && apt-get install -y  cmake libcairo2-dev libglpk-dev libgmp-dev libicu-dev libpng-dev libssl-dev libxml2-dev make pandoc zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN R -e 'remotes::install_cran("rlang")'
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("magrittr",upgrade="never", version = "2.0.3")'
RUN Rscript -e 'remotes::install_version("stringr",upgrade="never", version = "1.4.0")'
RUN Rscript -e 'remotes::install_version("tibble",upgrade="never", version = "3.1.7")'
RUN Rscript -e 'remotes::install_version("purrr",upgrade="never", version = "0.3.5")'
RUN Rscript -e 'remotes::install_version("dplyr",upgrade="never", version = "1.0.9")'
RUN Rscript -e 'remotes::install_version("colorspace",upgrade="never", version = "2.0-3")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never", version = "1.14.2")'
RUN Rscript -e 'remotes::install_version("ggplot2",upgrade="never", version = "3.3.6")'
RUN Rscript -e 'remotes::install_version("preprocessCore",upgrade="never", version = "1.58.0")'
RUN Rscript -e 'remotes::install_version("matrixStats",upgrade="never", version = "0.62.0")'
RUN Rscript -e 'remotes::install_version("hrbrthemes",upgrade="never", version = "0.8.0")'
RUN Rscript -e 'remotes::install_version("S4Vectors",upgrade="never", version = "0.34.0")'
RUN Rscript -e 'remotes::install_version("tidytable",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'remotes::install_version("edgeR",upgrade="never", version = "3.38.1")'
RUN Rscript -e 'remotes::install_version("statmod",upgrade="never", version = "1.4.36")'
RUN Rscript -e 'remotes::install_version("RANN",upgrade="never", version = "2.6.1")'
RUN Rscript -e 'remotes::install_version("SummarizedExperiment",upgrade="never", version = "1.26.1")'
RUN Rscript -e 'remotes::install_version("GGally",upgrade="never", version = "2.1.2")'
RUN Rscript -e 'remotes::install_version("lme4",upgrade="never", version = "1.1-29")'
RUN Rscript -e 'remotes::install_version("sva",upgrade="never", version = "3.44.0")'
RUN Rscript -e 'remotes::install_version("RUVnormalize",upgrade="never", version = "1.30.0")'
RUN Rscript -e 'remotes::install_version("RUVSeq",upgrade="never", version = "1.30.0")'
RUN Rscript -e 'remotes::install_version("WGCNA",upgrade="never", version = "1.71")'
RUN Rscript -e 'remotes::install_version("limma",upgrade="never", version = "3.52.1")'
RUN Rscript -e 'remotes::install_version("batchelor",upgrade="never", version = "1.12.0")'
RUN Rscript -e 'remotes::install_version("readr",upgrade="never", version = "2.1.2")'
RUN Rscript -e 'remotes::install_version("gtools",upgrade="never", version = "3.9.2.1")'
RUN Rscript -e 'remotes::install_version("igraph",upgrade="never", version = "1.3.1")'
RUN Rscript -e 'remotes::install_version("Biobase",upgrade="never", version = "2.56.0")'
RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN R -e 'remotes::install_local(upgrade="never")'
RUN rm -rf /build_zone
CMD R -e 'library(dockerfiler)'
