import numpy as np
from flask import Flask, jsonify, request
from werkzeug.exceptions import HTTPException
from flask_cors import CORS
from .constants import logger
import tensorflow as tf
import warnings
import os
import logging
import sys
from umap.parametric_umap import load_ParametricUMAP, ParametricUMAP
from sklearn import manifold, decomposition
from joblib import dump, load


# One of these should prevent tensorflow logs
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
warnings.filterwarnings("ignore", message=r"Passing", category=FutureWarning)
# logging.disable(logging.WARNING)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
logging.getLogger('tensorflow').setLevel(logging.FATAL)
tf.get_logger().setLevel(logging.ERROR)
# tf.autograph.set_verbosity(0)

# Create API
app = Flask(__name__)
# Enable CORS for the whole app
CORS(app)


@app.errorhandler(Exception)
def handle_exception(e):
    logger.exception(e)
    code = e.code if isinstance(e, HTTPException) else 500
    return jsonify({
        'message': str(e),
        'code': code
    }), code


def get_chembl_cddd():
    data = np.load('/_shared/chembl_cddd.npz', allow_pickle=True)
    smiles = data['smiles'].tolist()
    embedding = data['embedding']
    return smiles, embedding


_umap_model = None


def get_umap_model():
    global _umap_model
    if not _umap_model:
        logger.info('Loading UMAP model')
        # https://umap-learn.readthedocs.io/en/latest/parametric_umap.html
        try:
            # Load model
            _umap_model = load_ParametricUMAP('/_shared/p_umap')
            logger.info('Succesfully loaded UMAP model')
        except Exception as e:
            logger.exception('Loading UMAP model failed')
            raise e
            # sys.exit(1)
            # try:
            #     logger.exception('Retraining model')
            #     # Manually train the model
            #     chembl_full = np.load(
            #         '/home/data/chembl_ecfp4.npz', mmap_mode="r")['data']
            #     _umap_model = ParametricUMAP(verbose=True,
            #                             batch_size=512,
            #                             loss_report_frequency=1,
            #                             n_neighbors=200,
            #                             min_dist=0.7,
            #                             metric='dice')
            #     _umap_model.fit_transform(chembl_full[:10_000])
            #     _umap_model.save('/home/backend/models/p_umap')
            # except Exception:
            #     logger.exception('Training model failed, the API is not ready!')
    return _umap_model


# Singleton instance of the PCA model
_pca_model = None


def get_pca_model():
    global _pca_model
    if not _pca_model:
        try:
            # Try to load it from disk
            logger.info('Loading PCA model')
            _pca_model = load('/_shared/p_pca.joblib')
            logger.info('Succesfully loaded PCA model')
            return _pca_model
        except Exception:
            logger.exception('Loading UMAP model failed')
        # Otherwise fit it
        logger.info('Fitting PCA model')
        smiles, embedding = get_chembl_cddd()
        # Compute the model and projection, or only projection if the model is present
        if not _pca_model:
            # Create new model
            _pca_model = decomposition.PCA(n_components=2)
            # Fit the model on the embedding data
            _pca_model.fit(embedding)
            # Store it for later use
            dump(_pca_model, '/_shared/p_pca.joblib')
            logger.info('Succesfully fit and saved PCA model')
    return _pca_model


def get_method(method):
    if method == 'pca':
        model = get_pca_model()
    elif method == 'umap':
        model = get_umap_model()
    else:
        raise ValueError(f'Invalid projection method {method}')
    return model


@app.route("/api/projection/<method>/", methods=['GET'])
def get_umap(method):
    model = get_method(method)
    smiles, embedding = get_chembl_cddd()
    projection = model.transform(embedding).tolist()
    return jsonify(list({
        'structure': smiles[i],
        'projection': projection[i]
    } for i in range(len(smiles))))


@app.route("/api/similarity/", methods=['POST'])
def get_similarity():
    payload = request.get_json()
    smiles, embedding = get_chembl_cddd()
    reference = np.array(payload.get('reference'), dtype=np.float32)
    # Compute euclidian distance from ChEMBL to reference
    distance = ((embedding - reference)**2).sum(axis=1).tolist()
    return jsonify({ smiles[i]: distance[i] for i in range(len(smiles)) })


@app.route("/api/projection/<method>/", methods=['POST'])
def run_umap(method):
    model = get_method(method)
    payload = request.get_json()
    data = payload.get('data')
    if data is None or not isinstance(data, list) or len(data) == 0:
        raise ValueError('No valid data array in json body')

    data = np.array(data)
    return jsonify(model.transform(data).tolist())
