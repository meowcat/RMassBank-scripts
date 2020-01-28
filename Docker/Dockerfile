FROM bioconductor/devel_metabolomics2
MAINTAINER Michele Stravs stravsmi@eawag.ch

RUN mkdir -p /home/install
ADD setup.R /home/install/setup.R
RUN Rscript /home/install/setup.R

RUN cd /home/rstudio && git clone https://www.github.com/meowcat/RMassBank-scripts
