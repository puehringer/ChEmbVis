import umap
import numpy as np
from sklearn import manifold, decomposition
from ..constants import logger
from .tmap import tmap_projection, tmap_hash_projection
from .chembl_umap import umap_projection
from .chembl_tsne import tsne_projection


def compute_all_projections(data, additional):
    all_projections = ['pca', 'tsne', 'umap',
                       'chembl_umap', 'chembl_tsne', 'tmap']
    successful_projections = []
    projections = {}
    for t in all_projections:
        try:
            projections[t] = compute_projections(data, additional, {'type': t})
            successful_projections.append(t)
        except Exception:
            logger.exception(f'Error computing projection {t}')
    return successful_projections, projections


def compute_projections(data, additional, options):
    # Use embedding per default
    data_key = 'embedding'
    # TODO: Maybe create subclasses
    if options['type'] == 'pca':
        model = decomposition.PCA(n_components=2)
    elif options['type'] == 'tsne':
        model = manifold.TSNE(n_components=2, perplexity=30)
    elif options['type'] == 'mds':
        model = manifold.MDS(n_components=2, max_iter=100, n_init=1)
    elif options['type'] == 'umap':
        model = umap.UMAP()
    elif options['type'] == 'chembl_umap':
        def model(): return None
        model.fit_transform = umap_projection
        model.transform = umap_projection
    elif options['type'] == 'chembl_tsne':
        def model(): return None
        model.fit_transform = tsne_projection
        model.transform = tsne_projection
        data_key = 'smiles'
    elif options['type'] == 'tmap':
        def model(): return None
        model.fit_transform = tmap_hash_projection
        data_key = 'smiles'
    elif options['type'] == 'tmap_2':
        def model(): return None
        model.fit_transform = tmap_projection

    # Fit the model and project the data
    transformed_data = model.fit_transform(data[data_key])
    # Transform additional datapoints
    transformed_additional = model.transform(additional[data_key]) if additional and additional.get(
        data_key) and hasattr(model, 'transform') else None

    # Convert to list if we get an np array as output
    if type(transformed_data) is np.ndarray:
        transformed_data = transformed_data.tolist()
    if type(transformed_additional) is np.ndarray:
        transformed_additional = transformed_additional.tolist()
    return (transformed_data, transformed_additional, model)
