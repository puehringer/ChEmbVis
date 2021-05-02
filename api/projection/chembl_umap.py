import requests
import numpy as np
from typing import List
from ..constants import logger
from ..utils import cached
    
@cached
def get_projected_umap():
    try:
        return requests.get('http://api_umap:5000/api/projection/umap/').json()[:25_000]
    except Exception:
        logger.exception('Error fetching ChEMBL UMAP')


@cached
def get_projected_pca():
    try:
        return requests.get('http://api_umap:5000/api/projection/pca/').json()[:25_000]
    except Exception:
        logger.exception('Error fetching ChEMBL PCA')


@cached
def get_projected_pca_components():
    try:
        return np.array(requests.get('http://api_umap:5000/api/projection/pca/components').json())
    except Exception:
        logger.exception('Error fetching ChEMBL PCA Components')


def umap_projection(cddd_embedding: List[List[float]]):
    return requests.post('http://api_umap:5000/api/projection/umap/', json={'data': cddd_embedding}).json()


def pca_projection(cddd_embedding: List[List[float]]):
    return requests.post('http://api_umap:5000/api/projection/pca/', json={'data': cddd_embedding}).json()
