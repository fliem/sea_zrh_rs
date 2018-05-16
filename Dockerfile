FROM continuumio/miniconda3:latest
#RUN conda config --add channels conda-forge && \
RUN conda install -y pandas
RUN conda install -y joblib
RUN conda install -y seaborn
RUN conda install -y xlrd
RUN conda install scikit-learn
RUN pip install nilearn
RUN pip install -U feather-format

RUN apt-get update && apt-get install unzip
RUN mkdir -p /parcs && cd /parcs && \
    wget http://www.nil.wustl.edu/labs/petersen/Resources_files/Parcels.zip && \
    unzip Parcels.zip -d Gordon && rm -r Parcels.zip

RUN python -c "from nilearn import datasets;_=datasets.fetch_atlas_msdl()"
RUN python -c "from nilearn import datasets;_=datasets.fetch_atlas_yeo_2011()"

RUN mkdir -p /parcs/Schaefer
RUN cd /parcs/Schaefer && wget https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_17Networks_order.txt
RUN cd /parcs/Schaefer && wget https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_200Parcels_17Networks_order.txt
RUN cd /parcs/Schaefer && wget https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii.gz
RUN cd /parcs/Schaefer && wget https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm.nii.gz

RUN mkdir -p /code
COPY *.py /code/
RUN chmod +x /code/run.py

RUN python /code/yeo_sephem.py /parcs/Yeo_splithemi


ENTRYPOINT ["/code/run.py"]
