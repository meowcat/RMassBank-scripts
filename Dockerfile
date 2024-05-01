FROM bioconductor/bioconductor_docker

ADD . /home/rstudio/
RUN Rscript /home/rstudio/install.R

