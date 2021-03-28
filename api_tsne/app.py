from msp_tsne import ParametricTSNE
import numpy as np
from flask import Flask, jsonify, request
from werkzeug.exceptions import HTTPException
from flask_cors import CORS
from .constants import logger
import tensorflow as tf
import keras
import h5py
import warnings
import os
import logging
from tensorflow.python.keras.backend import set_session
import sys

# Create API
app = Flask(__name__)
# Enable CORS for the whole app
CORS(app)

try:
    # Why do I have to do this? https://github.com/tensorflow/tensorflow/issues/28287
    sess = tf.Session()
    graph = tf.get_default_graph()
    set_session(sess)
    # Load model
    p_tsne = ParametricTSNE(verbose=1, n_iter=200)
    p_tsne._build_model(2048, 2)
    p_tsne.model.load_weights('/_shared/p_tsne')
    p_tsne.model._make_predict_function()
    # graph = tf.get_default_graph()
    logger.info('Succesfully loaded TSNE model')
except Exception:
    logger.exception('Loading TNSE model failed')
    sys.exit(1)


@app.errorhandler(Exception)
def handle_exception(e):
    logger.exception(e)
    code = e.code if isinstance(e, HTTPException) else 500
    return jsonify({
        'message': str(e),
        'code': code
    }), code


@app.route("/api/tsne/", methods=['POST'])
def run_umap():
    payload = request.get_json()
    data = payload.get('data')
    if data is None or not isinstance(data, list) or len(data) == 0:
        raise ValueError('No valid data array in json body')

    with graph.as_default():
        set_session(sess)
        data = np.array(data).astype(np.float32)
        embedding = p_tsne.transform(data)
        return jsonify(embedding.tolist())
