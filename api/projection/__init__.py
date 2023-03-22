from typing import Any, Dict
import numpy as np
import pickle
import base64
from abc import ABC, abstractmethod
import requests
import json
from sklearn import decomposition, manifold
from ..constants import logger
from .tmap import tmap_projection, tmap_hash_projection
from .chembl_tsne import tsne_projection
from ..utils import catch_time
import math
import zlib
import umap

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




def all_models(data_key: str): 
    return {
        f'{data_key}_pca': lambda: SKLearnProjectionModel(data_key, decomposition.PCA(n_components=2)),
        f'{data_key}_tsne': lambda: SKLearnProjectionModel(data_key, manifold.TSNE(n_components=2, perplexity=30)),
        f'{data_key}_umap': lambda: SKLearnProjectionModel(data_key, umap.UMAP(n_components=2, n_neighbors=25, min_dist=0.5)),
        # f'{data_key}_umap': lambda: RemoteProjectionModel(data_key),
        # f'{data_key}_densmap': lambda: RemoteProjectionModel(data_key, model_options={'densmap': True}),
    }


models = {
    # **all_models('cddd'),
    # **all_models('graph'),
    **all_models('chemnet'),
    **all_models('molbert'),
    # **all_models('vae'),
    **all_models('maccs'),
    **all_models('ecfp0'),
    **all_models('ecfp2'),
    **all_models('ecfp4'),
    **all_models('ecfp6'),
    **all_models('rdkit'),
    **all_models('emb_clamp'),
    # **all_models('map4'),
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

    # Compute some predefined models for the precomputed embeddings
    precomputed_embeddings = options.get('precomputed_embeddings', [])
    for p in precomputed_embeddings:
        for type, model in [
            ('pca', SKLearnProjectionModel(p, decomposition.PCA(n_components=2))),
            # ('tsne', SKLearnProjectionModel(p, manifold.TSNE(n_components=2, perplexity=30))),
            ('umap', RemoteProjectionModel(p)),
            ('densmap', RemoteProjectionModel(p, model_options={'densmap': True})),
        ]:
            with catch_time(f'Computing another model for {p}'):
                try:
                    t = f'{p}_{type}'
                    projection = compute_projections(data, additional, {'type': t, 'model': model})
                    if projection[0] or projection[1]:
                        projections[t] = projection
                        successful_projections[t] = projection[3]
                except Exception:
                    logger.exception(f'Error computing another model for {p}')
    return successful_projections, projections


def trustworthiness_function():
    try:
        from sklearn.manifold import trustworthiness
        return trustworthiness
    except:
        pass
    try:
        from sklearn.manifold.t_sne import trustworthiness
        return trustworthiness
    except:
        pass
    logger.error('Error importing: from sklearn.manifold.t_sne import trustworthiness or from sklearn.manifold import trustworthiness')
    return None



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

    trustworthiness = trustworthiness_function()
    trustworthiness_score = trustworthiness(np.array(transformed_data), np.array(projected_data)) if trustworthiness and transformed_data else None
    trustworthiness_additional_score = trustworthiness(np.array(transformed_additional), np.array(projected_additional)) if trustworthiness and transformed_additional else None
    
    return (projected_data, projected_additional, model, {
        "model": model.serialize_model(),
        "trustworthiness": trustworthiness_score,
        "trustworthiness_additional": trustworthiness_additional_score,
        **model.information()
    })
