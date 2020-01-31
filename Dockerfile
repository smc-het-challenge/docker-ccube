FROM keyuan/docker-ccube
MAINTAINER "Kami Chiotti" chiotti@ohsu.edu

RUN install.r \
    mcclust

RUN installGithub.r \
    keyuan/ccube \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

COPY ./run_analysis_ccube_bb.R /home/pipeline/run_analysis_ccube_bb.R

RUN chmod +x /home/pipeline/run_analysis_ccube_bb.R

