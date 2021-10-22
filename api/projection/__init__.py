from typing import Any, Dict
import numpy as np
import pickle
import base64
from abc import ABC, abstractmethod
import requests
import json
from sklearn import decomposition, manifold
from sklearn.manifold.t_sne import trustworthiness
from ..constants import logger
from .tmap import tmap_projection, tmap_hash_projection
from .chembl_tsne import tsne_projection
from ..utils import catch_time
import math
import zlib

class ProjectionModel(ABC):
    # Default data key
    data_key = 'smiles'

    @classmethod
    def deserialize_model(cls, model_dump: str):
        return None

    def serialize_model(self) -> str:
        return None

    def information(self) -> Dict[str, Any]:
        return {}

    def distance_metric(self):
        return 'euclidean'

    @abstractmethod
    def fit_transform(self, data: any):
        pass

    def transform(self, data: any):
        return None

    def get_data(self, data: any = {}):
        if self.data_key not in data:
            raise ValueError(f'Invalid data_key {self.data_key}')
        return data.get(self.data_key)

class SKLearnProjectionModel(ProjectionModel):
    def __init__(self, data_key: str, model: any):
        super().__init__()
        self.data_key = data_key
        self.model = model

    def information(self) -> Dict[str, Any]:
        explained_variance = sum(self.model.explained_variance_ratio_[:2]) if hasattr(self.model, 'explained_variance_ratio_') else None
        return {
            "explained_variance": explained_variance if explained_variance is not None and not math.isnan(explained_variance) else None
        }

    @classmethod
    def deserialize_model(cls, model_dump: str):
        return SKLearnProjectionModel(**pickle.loads(zlib.decompress(base64.b64decode(model_dump))))

    def serialize_model(self) -> str:
        return base64.b64encode(zlib.compress(pickle.dumps({'model': self.model, 'data_key': self.data_key}), level=9)).decode('ascii')

    def fit_transform(self, data: any):
        return self.model.fit_transform(data)

    def transform(self, data: any):
        return self.model.transform(data) if getattr(self.model, 'transform') else None


class ChEMBLPCAModel(ProjectionModel):
    data_key = 'cddd'

    @classmethod
    def deserialize_model(cls, model_dump: str):
        return ChEMBLPCAModel()

    def serialize_model(self) -> str:
        return None

    def fit_transform(self, data: any):
        from .chembl_umap import pca_projection
        return pca_projection(data)

    def transform(self, data: any):
        from .chembl_umap import pca_projection
        return pca_projection(data)


class ChEMBLUMAPModel(ProjectionModel):
    data_key = 'cddd'

    @classmethod
    def deserialize_model(cls, model_dump: str):
        return ChEMBLUMAPModel()

    def serialize_model(self) -> str:
        return None

    def fit_transform(self, data: any):
        from .chembl_umap import umap_projection
        return umap_projection(data)

    def transform(self, data: any):
        from .chembl_umap import umap_projection
        return umap_projection(data)


class ChEMBLTSNEModel(ProjectionModel):

    @classmethod
    def deserialize_model(cls, model_dump: str):
        return ChEMBLTSNEModel()

    def serialize_model(self) -> str:
        return None

    def fit_transform(self, data: any):
        return tsne_projection(data.get('smiles'))

    def transform(self, data: any):
        return tsne_projection(data.get('smiles'))


class TMAPHashModel(ProjectionModel):

    def fit_transform(self, data: any):
        return tmap_hash_projection(data.get('smiles'))


class TMAPCDDDModel(ProjectionModel):
    data_key = 'cddd'

    def fit_transform(self, data: any):
        return tmap_projection(data)

class RemoteProjectionModel(ProjectionModel):
    def __init__(self, data_key: str, model_dump: str = None, model_options: any = None):
        super().__init__()
        self.data_key = data_key
        self.model_dump = model_dump
        self.model_options = model_options

    @classmethod
    def deserialize_model(cls, model_dump: str):
        return RemoteProjectionModel(**json.loads(model_dump))

    def serialize_model(self) -> str:
        return json.dumps({'model_dump': self.model_dump, 'data_key': self.data_key})


    def fit_transform(self, data: any):
        r = requests.get('http://api_umap:5000/api/projection/default_umap/', json={
            'data': data,
            'options': self.model_options
        })
        r.raise_for_status()
        result = r.json()
        self.model_dump = result['model_dump']
        return result['data']

    def transform(self, data: any):
        r = requests.get('http://api_umap:5000/api/projection/default_umap/', json={
            'data': data,
            'model_dump': self.model_dump,
            'options': self.model_options
        })
        r.raise_for_status()
        result = r.json()
        return result['data']


models = {
    # 'chembl_pca': lambda: ChEMBLPCAModel(),
    # 'chembl_tsne': lambda: ChEMBLTSNEModel(),
    # 'chembl_umap': lambda: ChEMBLUMAPModel(),
    'cddd_pca': lambda: SKLearnProjectionModel('cddd', decomposition.PCA(n_components=2)),
    # 'cddd_tsne': lambda: SKLearnProjectionModel('cddd', manifold.TSNE(n_components=2, perplexity=30)),
    'cddd_umap': lambda: RemoteProjectionModel('cddd'),
    'cddd_densmap': lambda: RemoteProjectionModel('cddd', model_options={'densmap': True}),
    # 'cddd_isomap': lambda: SKLearnProjectionModel('cddd', manifold.Isomap()),
    'graph_pca': lambda: SKLearnProjectionModel('graph', decomposition.PCA(n_components=2)),
    # 'graph_tsne': lambda: SKLearnProjectionModel('graph', manifold.TSNE(n_components=2, perplexity=30)),
    'graph_umap': lambda: RemoteProjectionModel('graph'),
    'graph_densmap': lambda: RemoteProjectionModel('graph', model_options={'densmap': True}),
    # 'graph_isomap': lambda: SKLearnProjectionModel('graph', manifold.Isomap()),
    # 'ecfp4_isomap': lambda: SKLearnProjectionModel('ecfp4', manifold.Isomap()),
    'chemnet_pca': lambda: SKLearnProjectionModel('chemnet', decomposition.PCA(n_components=2)),
    # 'chemnet_tsne': lambda: SKLearnProjectionModel('cddd', manifold.TSNE(n_components=2, perplexity=30)),
    'chemnet_umap': lambda: RemoteProjectionModel('chemnet'),
    'chemnet_densmap': lambda: RemoteProjectionModel('chemnet', model_options={'densmap': True}),
    'molbert_pca': lambda: SKLearnProjectionModel('molbert', decomposition.PCA(n_components=2)),
    # 'molbert_tsne': lambda: SKLearnProjectionModel('molbert', manifold.TSNE(n_components=2, perplexity=30)),
    'molbert_umap': lambda: RemoteProjectionModel('molbert'),
    'vae_pca': lambda: SKLearnProjectionModel('vae', decomposition.PCA(n_components=2)),
    # 'vae_tsne': lambda: SKLearnProjectionModel('vae', manifold.TSNE(n_components=2, perplexity=30)),
    'vae_umap': lambda: RemoteProjectionModel('vae'),
    'vae_densmap': lambda: RemoteProjectionModel('vae', model_options={'densmap': True}),
    'ecfp2_pca': lambda: SKLearnProjectionModel('ecfp2', decomposition.PCA(n_components=2)),
    # 'ecfp2_tsne': lambda: SKLearnProjectionModel('ecfp2', manifold.TSNE(n_components=2, perplexity=30)),
    'ecfp2_umap': lambda: RemoteProjectionModel('ecfp2'),
    'ecfp2_densmap': lambda: RemoteProjectionModel('ecfp2', model_options={'densmap': True}),
    'ecfp4_pca': lambda: SKLearnProjectionModel('ecfp4', decomposition.PCA(n_components=2)),
    # 'ecfp4_tsne': lambda: SKLearnProjectionModel('ecfp4', manifold.TSNE(n_components=2, perplexity=30)),
    'ecfp4_umap': lambda: RemoteProjectionModel('ecfp4'),
    'ecfp4_densmap': lambda: RemoteProjectionModel('ecfp4', model_options={'densmap': True}),
    'ecfp6_pca': lambda: SKLearnProjectionModel('ecfp6', decomposition.PCA(n_components=2)),
    # 'ecfp6_tsne': lambda: SKLearnProjectionModel('ecfp6', manifold.TSNE(n_components=2, perplexity=30)),
    'ecfp6_umap': lambda: RemoteProjectionModel('ecfp6'),
    'ecfp6_densmap': lambda: RemoteProjectionModel('ecfp6', model_options={'densmap': True}),
    # 'hash_tmap': lambda: TMAPHashModel(),
    # 'cddd_mds': lambda: SKLearnProjectionModel('cddd', manifold.MDS(n_components=2, max_iter=100, n_init=1)),
    # 'cddd_umap': lambda: SKLearnProjectionModel('cddd', umap.UMAP()),
    # 'tmap_2': lambda: TMAPCDDDModel(),
}


def compute_all_projections(data, additional, options={}):
    # data and additional is {'smiles': [...], 'cddd': [[...], [...], ...], 'molbert': [[...], [...], ...]}
    all_projections = options.get('all_projections', models.keys())
    successful_projections = {}
    projections = {}
    for t in all_projections:
        with catch_time(f'Computing {t}'):
            try:
                model = None
                if 'models' in options and t in options['models'] and t in models:
                    try:
                        model = models[t]().__class__.deserialize_model(
                            options['models'][t])
                    except:
                        logger.exception(
                            f'Model {t} could not be deserialized, skipping...')
                        continue

                projection = compute_projections(data, additional, {'type': t, 'model': model})
                if projection[0] or projection[1]:
                    projections[t] = projection
                    successful_projections[t] = projection[3]
            except Exception:
                logger.exception(f'Error computing projection {t}')
    return successful_projections, projections


def compute_projections(data, additional, options):
    projection_type = options['type']
    model = options.get('model')

    if not model:
        if projection_type in models:
            model = models.get(projection_type)()
        else:
            raise ValueError('Not implemented')

    transformed_data = model.get_data(data) if data else None
    transformed_additional = model.get_data(additional) if additional else None

    projected_data = model.fit_transform(transformed_data) if transformed_data else None
    projected_additional = model.transform(transformed_additional) if transformed_additional else None

    # Convert to list if we get an np array as output
    if type(projected_data) is np.ndarray:
        projected_data = projected_data.tolist()
    if type(projected_additional) is np.ndarray:
        projected_additional = projected_additional.tolist()

    trustworthiness_score = None
    try:
        trustworthiness_score = trustworthiness(np.array(transformed_data), np.array(projected_data)) if transformed_data else None
    except:
        logger.exception('Error computing trustworthiness')

    trustworthiness_additional_score = None
    try:
        trustworthiness_additional_score = trustworthiness(np.array(transformed_additional), np.array(projected_additional)) if transformed_additional else None
    except:
        logger.exception('Error computing additional trustworthiness')
    
    return (projected_data, projected_additional, model, {
        "model": model.serialize_model(),
        "trustworthiness": trustworthiness_score,
        "trustworthiness_additional": trustworthiness_additional_score,
        **model.information()
    })
