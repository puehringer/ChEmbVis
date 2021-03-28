import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
from flask import Flask, jsonify
from werkzeug.exceptions import HTTPException
from flask_cors import CORS
from flask_smorest import Api
from .constants import blp, logger
import warnings
import logging
import sys
# For some reason we need to import tensorflow at startup, as otherwise CDDD fails with ModuleNotFoundError: No module named 'tensorflow_core.keras'.
import tensorflow as tf
logger.info(f'GPU available: {tf.test.is_gpu_available()}')
from .api import InterpolationAPI  # noqa: F401
from .api import PSOAPI  # noqa: F401
from .api import ProjectionAPI  # noqa: F401
from .api import MoleculeImageAPI  # noqa: F401
from .api import EmbeddingAPI  # noqa: F401
from .api import MMPAPI  # noqa: F401

# Why do I have to do this?
sys.path.append('/home/backend')

# Create API
app = Flask(__name__)
# Enable CORS for the whole app
CORS(app)
# Specify the specification for flask_smorest
app.config['API_TITLE'] = 'API'
app.config['API_VERSION'] = 'v1'
app.config['OPENAPI_VERSION'] = '3.0.2'
app.config['OPENAPI_URL_PREFIX'] = '/api/spec'
app.config['OPENAPI_JSON_PATH'] = 'api.json'
app.config['OPENAPI_SWAGGER_UI_PATH'] = "/"
app.config['OPENAPI_SWAGGER_UI_URL'] = "https://cdn.jsdelivr.net/npm/swagger-ui-dist/"
# Create the flask_smorest api
api = Api(app)


@app.errorhandler(Exception)
def handle_exception(e):
    logger.exception(e)
    code = e.code if isinstance(e, HTTPException) else 500
    return jsonify({
        'message': str(e),
        'code': code
    }), code


api.register_blueprint(blp)
