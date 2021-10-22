FROM continuumio/miniconda3:4.9.2-alpine
ENV PATH /opt/conda/bin:$PATH

# Other deps
RUN conda install -q -y conda-forge::rdkit tensorflow
RUN conda install faerun
RUN conda install umap-learn
RUN conda install -c conda-forge notebook

# Fix tqdm
RUN conda install -c conda-forge ipywidgets
RUN jupyter nbextension enable --py widgetsnbextension

CMD ["jupyter", "notebook", "--notebook-dir=/home/jovyan/work", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]