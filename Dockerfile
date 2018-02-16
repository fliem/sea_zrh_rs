FROM continuumio/miniconda3:latest
RUN conda config --add channels conda-forge && \
    conda install -y pandas joblib nilearn
RUN conda install -y seaborn

RUN mkdir -p /code
COPY *.py /code/
RUN chmod +x /code/run.py

ENTRYPOINT ["/code/run.py"]
