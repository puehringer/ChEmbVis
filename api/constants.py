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
_cddd_model = None
def get_cddd_model():
    global _cddd_model
    if not _cddd_model:
        from cddd.inference import InferenceModel
        _cddd_model = InferenceModel('/_shared/p_cddd')
    return _cddd_model


# Create graph model for inference
_graph_model = None
def get_graph_model():
    global _graph_model
    if not _graph_model:
        # Load file and create lookup
        import pandas as pd
        df = pd.read_csv('/_shared/p_graph/muv_with_emb.csv')
        embeddings = df.loc[:, 'emb_0':'emb_255'].to_numpy()
        smiles_to_idx = {v: k for k, v in df.loc[:, 'smiles'].to_dict().items()}
        _graph_model = lambda smiles: embeddings[[smiles_to_idx[smi] for smi in smiles]]
    return _graph_model


# Create MolBERT model for inference
_molbert_model = None
def get_molbert_model():
    global _molbert_model
    if not _molbert_model:
        from molbert.utils.featurizer.molbert_featurizer import MolBertFeaturizer
        _molbert_model = MolBertFeaturizer('/_shared/p_molbert/checkpoints/last.ckpt', embedding_type='pooled', max_seq_len=200, device='cuda')
    return _molbert_model


# Create VAE model for inference
_vae_model = None
def get_vae_model():
    global _vae_model
    if not _vae_model:
        from chemvae.vae_utils import VAEUtils
        # from chemvae import mol_utils as mu
        _vae_model = VAEUtils(directory='/_shared/p_chemvae/zinc_properties')
    return _vae_model


