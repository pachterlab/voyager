FROM rocker/tidyverse:latest

ENV PROJ_VERSION="9.4.1"
ENV GDAL_VERSION="3.9.2"
ENV GEOS_VERSION="3.12.2"

COPY ./install_dev_osgeo.sh /rocker_scripts/install_dev_osgeo.sh
COPY ./install_geospatial.sh /rocker_scripts/install_geospatial.sh
RUN /rocker_scripts/install_dev_osgeo.sh
#COPY ./install_geospatial_unstable.sh /rocker_scripts/install_geospatial_unstable.sh
#RUN /rocker_scripts/install_geospatial_unstable.sh

# Bioc sysreqs
RUN apt-get update && apt-get install -y --no-install-recommends libfftw3-dev libmagick++-dev

RUN Rscript -e "BiocManager::install(version = 'devel', ask = FALSE)"
RUN Rscript -e "options(repos = c(CRAN='https://cloud.r-project.org')); BiocManager::install('Voyager', dependencies = TRUE)"
RUN Rscript -e "remotes::install_github('pachterlab/SpatialFeatureExperiment', ref='devel')"
RUN Rscript -e "remotes::install_github('pachterlab/voyager', ref='devel')"
