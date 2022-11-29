FROM python:3.9-slim-buster

WORKDIR /usr/src/app

COPY README.md pyproject.toml ./
ADD tsg ./tsg
ADD data ./data

RUN pip install .

ENTRYPOINT [ "transcript-generator" ]