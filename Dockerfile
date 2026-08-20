FROM continuumio/miniconda3:latest

WORKDIR /pipeline

COPY environment.yml .
RUN conda env update -n base -f environment.yml && conda clean -afy

COPY source/ source/