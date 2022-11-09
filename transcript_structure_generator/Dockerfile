FROM python:3-slim

WORKDIR /usr/src/app

COPY README.md requirements.txt pyproject.toml ./
ADD tsg ./tsg
ADD data ./data

RUN pip install .

CMD [ "transcript-generator", "--transcripts", "data/transcripts_short.csv", "--annotation", "data/Homo_sapiens.GRCh38.103.chr.gtf", "--prob_inclusion", "0.05" ]