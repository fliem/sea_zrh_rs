FROM continuumio/miniconda3:latest
RUN conda config --add channels conda-forge && \
    conda install -y pandas joblib nilearn
RUN conda install -y seaborn
RUN conda install -y xlrd

RUN apt-get update && apt-get install unzip
RUN mkdir -p /parcs && cd /parcs && \
    wget http://www.nil.wustl.edu/labs/petersen/Resources_files/Parcels.zip && \
    unzip Parcels.zip -d Gordon && rm -r Parcels.zip
RUN python -c "from nilearn import datasets;_=datasets.fetch_atlas_msdl()"

RUN mkdir -p /code
COPY *.py /code/
RUN chmod +x /code/run.py

ENTRYPOINT ["/code/run.py"]
