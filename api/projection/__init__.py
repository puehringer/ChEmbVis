import numpy as np
from sklearn import manifold, decomposition
from ..constants import logger
from .tmap import tmap_projection, tmap_hash_projection
from .chembl_umap import umap_projection, pca_projection
from .chembl_tsne import tsne_projection
from ..utils import catch_time, mol


def compute_all_projections(data, additional):
    all_projections = ['chembl_umap', 'chembl_pca', 'cddd_pca', 'hash_tmap', 'cddd_tsne', 'cddd_umap']
    successful_projections = []
    projections = {}
    for t in all_projections:
        with catch_time(f'Computing {t}'):
            try:
                projection = compute_projections(data, additional, {'type': t})
                if projection[0]:
                    projections[t] = projection
                    successful_projections.append(t)
            except Exception:
                logger.exception(f'Error computing projection {t}')
    return successful_projections, projections


def compute_projections(data, additional, options):
    # Use embedding per default
    data_key = 'embedding'
    # TODO: Maybe create subclasses
    if options['type'] == 'cddd_pca':
        model = decomposition.PCA(n_components=2)
    elif options['type'] == 'morgan_pca':
        model = decomposition.PCA(n_components=2)
        data_key = 'morgan'
    elif options['type'] == 'daylight_pca':
        model = decomposition.PCA(n_components=2)
        data_key = 'daylight'
    elif options['type'] == 'cddd_tsne':
        model = manifold.TSNE(n_components=2, perplexity=30)
    elif options['type'] == 'cddd_mds':
        model = manifold.MDS(n_components=2, max_iter=100, n_init=1)
    elif options['type'] == 'cddd_umap':
        # TODO: Importing on top causes weird memory/pointer issue.
        import umap
        model = umap.UMAP()
    elif options['type'] == 'chembl_umap':
        def model(): return None
        model.fit_transform = umap_projection
        model.transform = umap_projection
    elif options['type'] == 'chembl_pca':
        def model(): return None
        model.fit_transform = pca_projection
        model.transform = pca_projection
    elif options['type'] == 'chembl_tsne':
        def model(): return None
        model.fit_transform = tsne_projection
        model.transform = tsne_projection
        data_key = 'smiles'
    elif options['type'] == 'hash_tmap':
        def model(): return None
        model.fit_transform = tmap_hash_projection
        data_key = 'smiles'
    elif options['type'] == 'tmap_2':
        def model(): return None
        model.fit_transform = tmap_projection

    try:
        # Check if it is a valid fingerprint
        mol.to_fingerprint(data_key)
        # Fetch fingerprints if necessary
        data[data_key] = mol.mols_to_fingerprints(data['smiles'], data_key)
        if additional and additional.get('smiles'):
            additional[data_key] = mol.mols_to_fingerprints(additional['smiles'], data_key)
    except ValueError:
        # ignore
        pass

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
