# Chemical Embedding Exploration

This project is a (currently developed) application dedicated to the interactive analysis and exploration of chemical embeddings and projections. A live demo is available at http://cddd-explorer.pueh.xyz/. 

Several submodules make up the application:
* <a href="api/README.md"><b>api</b></a>: Main API exposing projections, <a target="_blank" href="https://github.com/jrwnter/mso">MSO</a>, rdkit, ...
* <a href="client/README.md"><b>client</b></a>: Interactive web application
* <a href="api_umap/README.md"><b>api_umap</b></a>: API for pretrained parametric UMAP
* <a href="api_tsne/README.md"><b>api_tsne</b></a>: API for pretrained parametric TSNE
* <a href="notebook/README.md"><b>notebook</b></a>: Jupyter notebook with several experiments

One reason why the API is split into several submodules is because they depend on different `tensorflow` versions, as the main API requires `tensorflow==1.15.0`, while the others depend on `tensorflow>=2.0`.

## Getting Started

The easiest way to run the whole application is by using docker[-compose]. To start local development, simply use `docker-compose start`. Please note that building the images might take a while, as we are using large dependencies such as `RDKit` and `Tensorflow`.

**GPU Support**: By default, the GPU enabled `tensorflow-gpu` is installed and docker-compose uses GPU capabilities for some containers (requires docker-compose 1.28+). For WSL2 users, please see https://docs.nvidia.com/cuda/wsl-user-guide/index.html first. For docker and docker-compose related issues, see https://docs.docker.com/config/containers/resource_constraints/#gpu and https://docs.docker.com/compose/gpu-support/#enabling-gpu-access-to-service-containers. 
If you want to use pure CPU implementations, please replace all `tensorflow-gpu` occurances with `tensorflow` in the corresponding `environment.yml` files. 
Also, you need to remove the `deploy` entries in the `docker-compose.yml`.

Please note that this project is actively developed such that no optimized builds/deployments are available yet.
