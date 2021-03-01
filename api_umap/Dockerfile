FROM continuumio/miniconda3:4.9.2-alpine
ENV PATH /opt/conda/bin:$PATH
ENV LANG C

# Copy everything into backend
WORKDIR /home/backend
COPY . ./
# Move the _shared folder to the image
RUN mv ./_shared /_shared

# Install conda dependencies
RUN conda update -n base conda
RUN conda env create -f environment.yml
RUN echo "source activate umap_api" > ~/.bashrc
ENV PATH /opt/conda/envs/umap_api/bin:$PATH
# Clean conda to reduce image size
RUN conda clean -tipsy

WORKDIR /home/backend/
CMD ["flask", "run", "--host=0.0.0.0"]