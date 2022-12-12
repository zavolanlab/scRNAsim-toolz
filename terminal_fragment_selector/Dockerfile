##### BASE #####
FROM python:3.10-slim-buster

##### VARIABLES #####
WORKDIR /Users/terminal-fragment-selector/

COPY requirements.txt /Users/terminal-fragment-selector/requirements.txt
COPY requirements_dev.txt /Users/terminal-fragment-selector/requirements_dev.txt

##### INSTALL #####
RUN apt-get update \
&& apt-get install gcc -y \
&& apt-get clean
RUN pip install -r requirements.txt -r requirements_dev.txt




