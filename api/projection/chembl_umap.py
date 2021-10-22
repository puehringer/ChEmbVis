import requests
import numpy as np
from typing import List
from ..constants import logger
from ..utils import cached


def _json_or_raise_status(req):
    req.raise_for_status()
    return req.json()

    
@cached
def get_projected_umap():
    try:
        return _json_or_raise_status(requests.get('http://api_umap:5000/api/projection/umap/'))[:10_000]
    except Exception:
        logger.exception('Error fetching ChEMBL UMAP')


@cached
def get_projected_pca():
    try:
        return _json_or_raise_status(requests.get('http://api_umap:5000/api/projection/pca/'))[:10_000]
    except Exception:
        logger.exception('Error fetching ChEMBL PCA')


@cached
def get_projected_pca_components():
    try:
        return np.array(_json_or_raise_status(requests.get('http://api_umap:5000/api/projection/pca/components')))
    except Exception:
        logger.exception('Error fetching ChEMBL PCA Components')


def umap_projection(cddd_embedding: List[List[float]]):
    return _json_or_raise_status(requests.post('http://api_umap:5000/api/projection/umap/', json={'data': cddd_embedding}))


def pca_projection(cddd_embedding: List[List[float]]):
    return _json_or_raise_status(requests.post('http://api_umap:5000/api/projection/pca/', json={'data': cddd_embedding}))
