FROM python:3.11.1-slim-bullseye

# MAINTAINER add_maintainer

# set names for user, group and user home
ARG USER="bioc"
ARG GROUP="bioc"
ARG WORKDIR="/home/${USER}"

# create user, group, user home & makes the user and group
# the owner of the user home
RUN mkdir -p $WORKDIR \
  && groupadd -r $GROUP \
  && useradd --no-log-init -r -g $GROUP $USER \
  && chown -R ${USER}:${GROUP} $WORKDIR \
  && chmod 700 $WORKDIR

# set the user, make sure the location where pip
# installs binaries/executables is part of the $PATH
# and set the working directory to the user's home
USER $USER
ENV PATH="${WORKDIR}/.local/bin:${PATH}"
WORKDIR $WORKDIR

# copy entire content of the current directory to the
# working directory; ensure that this does not contain any
# files that you don't want to have end up in the Docker
# image, especially not any secrets!
# use `.dockerignore` to set files that you do NOT want
# Docker to copy
COPY --chown=${USER}:${GROUP} . $WORKDIR

# install app and development dependencies
# assumes that dependencies in `requirements.txt` are
# automatically installed when installing the app; if
# that is not the case, they need to be installed
# _before_ installing the app
RUN pip install . \
  && pip install --no-cache-dir -r requirements-dev.txt

# set default command; optional
CMD ["readsequencer"]
