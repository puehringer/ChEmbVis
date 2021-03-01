import requests
from typing import List
from ..constants import logger

# Fetch the precomputed ChEMBL UMAP and PCA on start
projected_umap = None
try:
    projected_umap = requests.get('http://api_umap:5000/api/umap/').json()[:25_000]
except Exception:
    logger.exception('Error fetching ChEMBL UMAP')

projected_pca = None
try:
    projected_pca = requests.get('http://api_umap:5000/api/pca/').json()[:25_000]
except Exception:
    logger.exception('Error fetching ChEMBL PCA')


def umap_projection(cddd_embedding: List[List[float]]):
    return requests.post('http://api_umap:5000/api/umap/', json={'data': cddd_embedding}).json()


def pca_projection(cddd_embedding: List[List[float]]):
    return requests.post('http://api_umap:5000/api/pca/', json={'data': cddd_embedding}).json()
