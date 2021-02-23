from ..utils.mol import _string_to_mol, AllChem
import requests
import numpy as np


def tsne_projection(mols):
    # Convert to mols
    all_mols = {s: _string_to_mol(s) for s in mols}
    # Compute fingerprints where possible
    fingerprints = {s: np.array(AllChem.GetMorganFingerprintAsBitVect(
        m, radius=2, nBits=2048)).tolist() for s, m in all_mols.items() if m is not None}
    # Use API to compute embedding
    embedding = requests.post('http://api_tsne:5000/api/tsne/',
                              json={'data': list(fingerprints.values())}).json()
    # Map embeddings back to SMILES
    lookup = {s: embedding[i] for i, s in enumerate(fingerprints.keys())}
    # Return w.r.t. all mols
    return [lookup.get(s) for s in mols]
