import logging
from flask_smorest import Blueprint

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s")
logger = logging.getLogger('api')
logger.setLevel(logging.INFO)

# Blueprint for the API
blp = Blueprint(
    'api', 'api', url_prefix='/api'
)

# Create CDDD model for inference
inference_model = None
def get_inference_model():
    global inference_model
    if not inference_model:
        from cddd.inference import InferenceModel
        inference_model = InferenceModel('/_shared/p_cddd')
    return inference_model
