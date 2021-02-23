FROM continuumio/miniconda3:4.9.2
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
RUN echo "source activate api" > ~/.bashrc
ENV PATH /opt/conda/envs/api/bin:$PATH
# Clean conda to reduce image size
RUN conda clean -tipsy

# Install CDDD
WORKDIR /home/lib/
RUN git clone https://github.com/jrwnter/cddd.git
RUN pip install cddd/
# WORKDIR /home/lib/cddd
# RUN chmod +x download_default_model.sh
# RUN ./download_default_model.sh
# TODO: Some hack here
# RUN cp default_model /opt/conda/envs/api/lib/python3.7/site-packages/cddd/data/

# Install MSO
WORKDIR /home/lib/
RUN git clone https://github.com/jrwnter/mso.git
RUN pip install mso/

# Add to path
ENV PATH /home/lib:$PATH

WORKDIR /home/backend/
CMD ["flask", "run", "--host=0.0.0.0", "--eager-loading"]
