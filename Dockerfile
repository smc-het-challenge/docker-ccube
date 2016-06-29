FROM rocker/hadleyverse
MAINTAINER "Ke Yuan" ke.yuan.09@gmail.com

RUN apt-get -qq update \
		&& apt-get install --no-install-recommends -y \
			libcurl4-openssl-dev \
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


RUN install2.r --error \
	doParallel \
	foreach

RUN installGithub.r \
    keyuan/ccube \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN mkdir /home/pipeline

COPY ./create_ccfclust_inputs.py /home/pipeline/create_ccfclust_inputs.py
COPY ./run_analysis_ccube.R /home/pipeline/run_analysis_ccube.R
COPY ./run_purity.R /home/pipeline/run_purity.R

RUN chmod +x /home/pipeline/create_ccfclust_inputs.py \
    && chmod +x /home/pipeline/run_analysis_ccube.R \
    && chmod +x /home/pipeline/run_purity.R



ENV PATH=/home/pipeline:$PATH

## ENTRYPOINT ["Rscript", "/home/run_analysis_ccube.R"]
