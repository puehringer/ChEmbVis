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


# https://umap-learn.readthedocs.io/en/latest/parametric_umap.html
try:
    # Load model
    p_umap = load_ParametricUMAP('/_shared/p_umap')
    logger.info('Succesfully loaded UMAP model')
except Exception:
    logger.exception('Loading UMAP model failed')
    sys.exit(1)
    # try:
    #     logger.exception('Retraining model')
    #     # Manually train the model
    #     chembl_full = np.load(
    #         '/home/data/chembl_ecfp4.npz', mmap_mode="r")['data']
    #     p_umap = ParametricUMAP(verbose=True,
    #                             batch_size=512,
    #                             loss_report_frequency=1,
    #                             n_neighbors=200,
    #                             min_dist=0.7,
    #                             metric='dice')
    #     p_umap.fit_transform(chembl_full[:10_000])
    #     p_umap.save('/home/backend/models/p_umap')
    # except Exception:
    #     logger.exception('Training model failed, the API is not ready!')


@app.errorhandler(Exception)
def handle_exception(e):
    logger.exception(e)
    code = e.code if isinstance(e, HTTPException) else 500
    return jsonify({
        'message': str(e),
        'code': code
    }), code


@app.route("/api/umap/", methods=['GET'])
def get_umap():
    data = np.load('/_shared/chembl_cddd.npz', allow_pickle=True)
    smiles = data['smiles'].tolist()
    embedding = data['embedding']
    projection = p_umap.transform(embedding).tolist()
    return jsonify(list({
        'structure': smiles[i],
        'projection': projection[i]
    } for i in range(len(smiles))))


@app.route("/api/umap/", methods=['POST'])
def run_umap():
    payload = request.get_json()
    data = payload.get('data')
    if data is None or not isinstance(data, list) or len(data) == 0:
        raise ValueError('No valid data array in json body')

    data = np.array(data)
    return jsonify(p_umap.transform(data).tolist())
