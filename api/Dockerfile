FROM continuumio/miniconda3:4.9.2
ENV PATH /opt/conda/bin:$PATH
ENV LANG C

WORKDIR /home/backend
# Copy environment into backend
COPY ./environment.yml ./environment.yml

# Install conda dependencies
RUN conda update -n base conda
RUN conda env create -f environment.yml
RUN echo "source activate api" > ~/.bashrc
ENV PATH /opt/conda/envs/api/bin:$PATH
# Clean conda to reduce image size
RUN conda clean -tipsy

# Copy everything else
COPY . ./
# Move the _shared folder to the image
RUN mv ./_shared /_shared

# HACK: For some reason chemical_vae installs Keras 2.0.7, which causes https://github.com/keras-team/keras/issues/7736 
RUN sed -i 's/v.constraint = constraint//' /opt/conda/envs/api/lib/python3.7/site-packages/keras/backend/tensorflow_backend.py

# Optionally configure Jupyter Notebook
RUN jupyter nbextension enable --py widgetsnbextension

WORKDIR /home/backend/
CMD ["flask", "run", "--host=0.0.0.0", "--eager-loading"]
