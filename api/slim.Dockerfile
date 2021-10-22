FROM python:3.7-buster

WORKDIR /home/backend

RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && apt-get purge -y --auto-remove \
    && rm -rf /var/lib/apt/lists/*

# Copy environment into backend
COPY ./requirements.txt ./requirements.txt

# Install latest pip and setuptools + cython to allow building of hdbscan
RUN python -m pip install --upgrade pip setuptools cython

RUN pip install -r ./requirements.txt --no-cache-dir
# Install hdbscan without binaries, otherwise a numpy binary issue arises: https://github.com/scikit-learn-contrib/hdbscan/issues/457
RUN pip install hdbscan --no-cache-dir --no-binary :all: --no-build-isolation

# Copy everything else
COPY . ./
# Move the _shared folder to the root of the image
RUN mv ./_shared /_shared

# HACK: For some reason chemical_vae installs Keras 2.0.7, which causes https://github.com/keras-team/keras/issues/7736 
# RUN sed -i 's/v.constraint = constraint//' /opt/conda/envs/api/lib/python3.7/site-packages/keras/backend/tensorflow_backend.py

WORKDIR /home/backend/
CMD ["flask", "run", "--host=0.0.0.0", "--eager-loading"]