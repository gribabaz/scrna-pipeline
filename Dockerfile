FROM continuumio/miniconda3:latest

WORKDIR /pipeline

COPY environment.yml .
RUN conda env update -n base -f environment.yml && conda clean -afy

ENV PATH=/opt/conda/envs/scrna-pipeline/bin:$PATH

COPY source/ source/