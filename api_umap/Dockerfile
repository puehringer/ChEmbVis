FROM continuumio/miniconda3:4.9.2-alpine
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

WORKDIR /home/backend/
CMD ["flask", "run", "--host=0.0.0.0", "--eager-loading"]