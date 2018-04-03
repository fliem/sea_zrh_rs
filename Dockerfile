FROM continuumio/miniconda3:latest
#RUN conda config --add channels conda-forge && \
RUN conda install -y pandas
RUN conda install -y joblib
RUN conda install -y seaborn
RUN conda install -y xlrd
RUN conda install scikit-learn
RUN pip install nilearn

RUN apt-get update && apt-get install unzip
RUN mkdir -p /parcs && cd /parcs && \
    wget http://www.nil.wustl.edu/labs/petersen/Resources_files/Parcels.zip && \
    unzip Parcels.zip -d Gordon && rm -r Parcels.zip
RUN python -c "from nilearn import datasets;_=datasets.fetch_atlas_msdl()"
RUN python -c "from nilearn import datasets;_=datasets.fetch_atlas_basc_multiscale_2015(version='sym')"

RUN mkdir -p /code
COPY *.py /code/
RUN chmod +x /code/run.py

ENTRYPOINT ["/code/run.py"]
