FROM continuumio/miniconda3

WORKDIR /usr/src/app

COPY . ./
RUN --mount=type=cache,target=/opt/conda/pkgs conda install -c anaconda --file requirements.txt
CMD [ "python", "main-fast.py"]