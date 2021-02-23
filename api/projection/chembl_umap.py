import requests
from typing import List
from ..constants import logger

# Fetch the precomputed ChEMBL UMAP on start
projected_umap = None
try:
    projected_umap = requests.get('http://api_umap:5000/api/umap/').json()
except Exception:
    logger.exception('Error fetching ChEMBL UMAP')


def umap_projection(cddd_embedding: List[List[float]]):
    # Convert to mols
    # all_mols = {s: _string_to_mol(s) for s in mols}
    # Compute fingerprints where possible
    # fingerprints = {s: np.array(AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=512)).tolist() for s, m in all_mols.items() if m is not None}
    # Use API to compute embedding
    # embedding = requests.post('http://api_umap:5000/api/umap/', json={'data': list(fingerprints.values())}).json()
    # Map embeddings back to SMILES
    # lookup = {s: embedding[i] for i, s in enumerate(fingerprints.keys())}
    # Return w.r.t. all mols
    # return [lookup.get(s) for s in mols]
    return requests.post('http://api_umap:5000/api/umap/', json={'data': cddd_embedding}).json()
