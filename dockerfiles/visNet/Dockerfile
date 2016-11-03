FROM ubuntu:14.04

RUN sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list' \
    && gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 \
    && gpg -a --export E084DAB9 | sudo apt-key add -

# [1] Get noninteractive frontend for Debian to avoid some problems:
#    debconf: unable to initialize frontend: Dialog
# [2] Always combine RUN apt-get update with apt-get install in the same RUN statement for cache busting
# https://docs.docker.com/engine/articles/dockerfile_best-practices/

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections \
    && apt-get update -y \
    && apt-get install -y \
        gcc \
        g++ \
        zlib1g-dev \
        libcurl4-openssl-dev \
        r-base \ 
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 


RUN echo 'source("http://bioconductor.org/biocLite.R"); biocLite(c("supraHex", "graph", "Rgraphviz", "Biobase"));q()'|R --no-save

RUN echo 'install.packages("devtools", repos="http://cran.r-project.org", clean=TRUE);q()' | R --no-save \
    && echo 'devtools::install_version("data.table", "1.9.6", repos="http://cran.r-project.org", clean=TRUE);q()' | R --no-save \
    && echo 'devtools::install_version("igraph", "1.0.1", repos="http://cran.r-project.org", clean=TRUE);q()' | R --no-save \
    && echo 'devtools::install_version("dnet", "1.0.9", repos="http://cran.r-project.org", clean=TRUE);q()' | R --no-save \
    && echo 'devtools::install_version("argparse", "1.0.1", repos="http://cran.r-project.org", clean=TRUE);q()' | R --no-save \
    && echo 'devtools::install_version("jsonlite", "1.0", repos="http://cran.r-project.org", clean=TRUE);q()' | R --no-save

COPY Rscripts/* /scripts/
COPY Inputs/* /scripts/  