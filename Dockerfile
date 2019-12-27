FROM ubuntu
MAINTAINER "Ke Yuan" ke.yuan.09@gmail.com

ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get -qq update \
        && apt-get install --no-install-recommends -y \
            libcurl4-openssl-dev \ 
            libssl-dev \
            apt-transport-https \
            python-dev \
            libc-dev \
            python-pip \
            pkg-config \
            liblzma-dev \
            libbz2-dev \
            libpcre3-dev \
            build-essential \
            libblas-dev \
            liblapack-dev \
            gfortran \
            libzmq3-dev \
            curl \
            libfreetype6-dev \
            libpng-dev \
            net-tools \
            procps \
            libreadline-dev \
        && pip install distribute --upgrade \
        && pip install \
            cython \
            pysam \ 
            pyvcf \
        && apt-get autoremove -y \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get -qq update \
    && apt-get install -y --no-install-recommends \
        littler \
        r-base \
        r-base-dev \
        r-recommended \
        && echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
    && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("xml2")' \
    && Rscript -e 'install.packages("rversions")' \
    && Rscript -e 'install.packages("roxygen2")' \
    && Rscript -e 'install.packages("dplyr")' \
    && Rscript -e 'install.packages("mcclust")' \
    && Rscript -e 'install.packages("doParallel")' \
    && Rscript -e 'install.packages("foreach")'\
    && Rscript -e 'install.packages("devtools")'

RUN Rscript -e 'devtools::install_github("keyuan/ccube")' \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN install.r \
    colorspace \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN mkdir /home/pipeline

COPY ./create_ccfclust_inputs.py /home/pipeline/create_ccfclust_inputs.py
COPY ./run_analysis_ccube.R /home/pipeline/run_analysis_ccube_test.R
COPY ./run_analysis_ccube.R /home/pipeline/run_analysis_ccube.R
COPY ./run_analysis_ccube.R /home/pipeline/run_analysis_ccube_bb.R


RUN chmod +x /home/pipeline/create_ccfclust_inputs.py \
    && chmod +x /home/pipeline/run_analysis_ccube_test.R \
    && chmod +x /home/pipeline/run_analysis_ccube.R \
    && chmod +x /home/pipeline/run_analysis_ccube_bb.R   
  
ENV PATH=/home/pipeline:$PATH
