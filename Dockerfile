FROM continuumio/miniconda3

WORKDIR /usr/src/app

COPY . ./
RUN --mount=type=cache,target=/opt/conda/pkgs conda install -c anaconda --file requirements.txt
ENV FLASK_DEBUG=1
CMD [ "flask", "run", "--port", "5000", "--host", "0.0.0.0"]