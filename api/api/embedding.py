from flask import jsonify
from flask_smorest import abort
import sklearn.neighbors as neighbors
from sklearn.metrics import pairwise_distances
import numpy as np
from sklearn.utils import graph
import hdbscan
import joblib
from ..schema import EmbeddingArgsSchema, ServerCollectionSchema
from ..utils import cached, mol, catch_time, parallelized
from ..constants import blp, logger
from ..projection import compute_all_projections
from ..projection.chembl_umap import get_projected_umap, get_projected_pca
from ..models import embedding_models

@blp.route('/embedding/', methods=['GET'])
@cached
def get_embedding():
    # TODO: Move to different API
    if not get_projected_umap() or not get_projected_pca():
        abort(404)
    structures = list(
        map(lambda p: p.get('structure'), get_projected_umap()))
    umap_projections = list(
        map(lambda p: p.get('projection'), get_projected_umap()))
    pca_projections = list(
        map(lambda p: p.get('projection'), get_projected_pca()))
    properties = list(parallelized(mol.compute_properties, structures))

    return jsonify([{
        'structure': structure,
        'properties': properties[i],
        'projection': {
            'chembl_umap': umap_projections[i],
            'chembl_pca': pca_projections[i]
        },
    } for i, structure in enumerate(structures)])


@blp.route('/embedding/', methods=['POST'])
@blp.arguments(EmbeddingArgsSchema, location='json')
@blp.response(200, ServerCollectionSchema)
@cached
def compute_embedding(args):
    raw_structures = args.get('structures')
    include_embedding = args.get('include_embedding')

    with joblib.parallel_backend('loky', n_jobs=2):

        with catch_time(f'Filtering {len(raw_structures)} structures'):
            # structures = [s for s in structures if mol.is_valid_mol(s)]
            # Keep a preprocessed -> original map to send back the original structure too
            structure_map = {new: original for (original, new) in map(lambda s: (s, mol.preprocess_smiles(s['smiles'])), raw_structures) if mol.is_valid_mol(new)}
            structures = list(structure_map.keys())

        computed_embeddings = {}

        nearest_neighbors = [{} for _ in range(len(structures))]
        def compute_nearest_neighbors(data, metric, key: str):
            neigh = neighbors.NearestNeighbors(n_neighbors=min(50, len(data) - 1), metric=metric)
            neigh.fit(data)
            all_dist, all_ind = neigh.kneighbors(return_distance=True)
            all_dist = all_dist.tolist()
            all_ind = all_ind.tolist()

            for i in range(len(nearest_neighbors)):
                nearest_neighbors[i][key] = {
                    'distance_metric': 'function' if callable(metric) else metric,
                    'knn_dist': all_dist[i],
                    'knn_ind': all_ind[i]
                }
                # nearest_neighbors[i][f'{key}_geodesic_5'] = {
                #     'distance_metric': 'function' if callable(distance_metric) else distance_metric,
                #     'knn_dist': dist_matrix_5[i].tolist()[1:51],
                #     'knn_ind': np.argsort(dist_matrix_5[i]).tolist()[1:51]
                # }
                # nearest_neighbors[i][f'{key}_geodesic_15'] = {
                #     'distance_metric': 'function' if callable(distance_metric) else distance_metric,
                #     'knn_dist': dist_matrix_15[i].tolist()[1:51],
                #     'knn_ind': np.argsort(dist_matrix_15[i]).tolist()[1:51]
                # }
                # nearest_neighbors[i][f'{key}_geodesic_50'] = {
                #     'distance_metric': 'function' if callable(distance_metric) else distance_metric,
                #     'knn_dist': dist_matrix_50[i].tolist()[1:51],
                #     'knn_ind': np.argsort(dist_matrix_50[i]).tolist()[1:51]
                # }

            # kng_5 = neighbors.kneighbors_graph(neigh, n_neighbors=min(5, len(data) - 1), metric=distance_metric, mode='distance')
            # dist_matrix_5 = graph.graph_shortest_path(kng_5, method='auto', directed=False)
            # kng_15 = neighbors.kneighbors_graph(neigh, n_neighbors=min(15, len(data) - 1), metric=distance_metric, mode='distance')
            # dist_matrix_15 = graph.graph_shortest_path(kng_15, method='auto', directed=False)
            # kng_50 = neighbors.kneighbors_graph(neigh, n_neighbors=min(50, len(data) - 1), metric=distance_metric, mode='distance')
            # dist_matrix_50 = graph.graph_shortest_path(kng_50, method='auto', directed=False)


        clusters = [{} for _ in range(len(structures))]
        def compute_clusters(data, metric, key: str):
            clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size=20, min_samples=10)
            clusterer.fit(data)
            labels = clusterer.labels_.tolist()
            
            for i in range(len(clusters)):
                clusters[i][key] = {
                    'distance_metric': 'function' if callable(metric) else metric,
                    'label': str(labels[i])
                }
        # Compute everything for precomputed embeddings
        # TODO: Combine with below function (as it is duplicated) 
        try:
            precomputed_embeddings = (list(structure_map.values())[0].get('embeddings') or {}).keys()
            logger.info(f'Processing existing embeddings: {precomputed_embeddings}')
            for key in precomputed_embeddings:
                with catch_time(f'Computing {key} embedding for {len(structures)} structures'):
                    computed_embedding = np.array([s['embeddings'][key] for s in structure_map.values()])
                with catch_time(f'Computing pairwise distances for {key} (single job)'):
                    distances = pairwise_distances(computed_embedding, metric='euclidean', n_jobs=1)
                # with catch_time(f'Computing pairwise distances for {key} (multiple jobs)'):
                #     distances = pairwise_distances(computed_embedding, metric=embedding.distance_metric)
                with catch_time(f'Computing knn for {key}'):
                    compute_nearest_neighbors(distances, metric='precomputed', key=key)
                with catch_time(f'Computing clusters for {key}'):
                    compute_clusters(distances.astype(np.float64), metric='precomputed', key=key)

                computed_embeddings[key] = computed_embedding.tolist()
        except Exception:
            logger.exception(f'Error computing properties of passed embeddings')

        for key, embedding in embedding_models.items():
            try:
                if embedding.can_encode:
                    with catch_time(f'Computing {key} embedding for {len(structures)} structures'):
                        computed_embedding = embedding.encode(structures if not embedding.use_raw_smiles else [structure_map[smi]['smiles'] for smi in structures])
                    with catch_time(f'Computing pairwise distances for {key} (single job)'):
                        distances = pairwise_distances(computed_embedding, metric=embedding.distance_metric, n_jobs=1)
                    # with catch_time(f'Computing pairwise distances for {key} (multiple jobs)'):
                    #     distances = pairwise_distances(computed_embedding, metric=embedding.distance_metric)
                    with catch_time(f'Computing knn for {key}'):
                        compute_nearest_neighbors(distances, metric='precomputed', key=key)
                    with catch_time(f'Computing clusters for {key}'):
                        compute_clusters(distances.astype(np.float64), metric='precomputed', key=key)

                    computed_embeddings[key] = computed_embedding.tolist()
            except Exception:
                logger.exception(f'Error computing embedding for {key}')

        with catch_time(f'Computing projections for {len(structures)} structures'):
            projections, projection = compute_all_projections({'smiles': structures, **computed_embeddings }, None, {
                'precomputed_embeddings': precomputed_embeddings
            })

        # for key, value in projection.items():
        #     with catch_time(f'Computing pairwise distances for {key}'):
        #         distances = pairwise_distances(value[0], metric='euclidean')
        #     with catch_time(f'Computing knn for {key}'):
        #         compute_nearest_neighbors(distances, metric='precomputed', key=key)
        #     with catch_time(f'Computing clusters for {key}'):
        #         compute_clusters(distances.astype(np.float64), metric='precomputed', key=key)

    result = {
        'data': [{
            'structure': structure,
            'original_structure': structure_map[structure]['smiles'],
            # 'decoded': decoded_smiles_list[i],
            'embedding': { key: embedding[i] for key, embedding in computed_embeddings.items() } if include_embedding else None,
            'projection': {t: projection[t][0][i] for t in projections},
            'nearest_neighbors': nearest_neighbors[i],
            'clusters': clusters[i],
            'properties': {**mol.compute_properties(structure)} # , "egfr": egfr_scores[i], "bace": bace_scores[i]
        } for i, structure in enumerate(structures)],
        "projections": {} # projections Do not include projections due to size limits
    }

    from datetime import datetime
    import json
    with open(f'/_shared/{datetime.now().strftime("%Y%m%d_%H%M%S")}.json', 'w') as file:
        logger.info(f'Storing file as {file.name}')
        json.dump(result, file)

    return result
