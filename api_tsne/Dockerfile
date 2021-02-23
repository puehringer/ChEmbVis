FROM continuumio/miniconda3:4.9.2-alpine
ENV PATH /opt/conda/bin:$PATH

# Copy everything into backend
WORKDIR /home/backend
COPY . ./

# Install conda dependencies
RUN conda update -n base conda
RUN conda env create -f environment.yml
RUN echo "source activate api_tsne" > ~/.bashrc
ENV PATH /opt/conda/envs/api_tsne/bin:$PATH

EXPOSE 5000

WORKDIR /home/backend/
CMD ["flask", "run", "--host=0.0.0.0"]