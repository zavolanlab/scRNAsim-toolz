###### BASE IMAGE ######
FROM continuumio/miniconda3:4.12.0

####### METADATA #######
LABEL base_image="continuumio/miniconda3:4.12.0"
LABEL version="1.0"
LABEL software="scRNAsim-toolz"
LABEL software.version="v0.1.1"
LABEL about.home="https://github.com/zavolanlab/scRNAsim-toolz"

###### MAINTAINER ######
LABEL maintainer="Máté Balajti <mate.balajti@unibas.ch>"

##### INSTALLATION #####

# COPY THE YAML & INSTALL SOFTWARE WITH CONDA
WORKDIR /usr/src/app
COPY ./ ./
RUN conda env create --file environment.yml \
    && conda clean --all

# VARIABLES
ARG WORKDIR="/home/USER"
ARG USER="USER"
ARG GROUP="GROUP"
ENV PATH="${WORKDIR}:${PATH}"

# CREATE USER
RUN groupadd -r ${GROUP} && useradd --no-log-init -r -g ${GROUP} ${USER}

# SET ENVIRONMENT
WORKDIR ${WORKDIR}
RUN chown -R ${USER}:${GROUP} ${WORKDIR} && chmod 700 ${WORKDIR}
USER ${USER}
RUN echo "source activate scrnasim-toolz" > ~/.bashrc
ENV PATH /opt/conda/envs/scrnasim-toolz/bin:$PATH

# SET ENTRYPOINT
CMD ["transcript-sampler", "-v"]