import numpy as np
import tmap
from mhfp.encoder import MHFPEncoder
from ..utils.mol import _string_to_mol
from ..utils import logger


def tmap_hash_projection(mols):
    all_mols = [_string_to_mol(mol) for mol in mols]
    mols = list(filter(None, all_mols))

    # The number of permutations used by the MinHashing algorithm
    perm = 512

    # Initializing the MHFP encoder with 512 permutations
    enc = MHFPEncoder(perm)

    # Create MHFP fingerprints from mol
    # The fingerprint vectors have to be of the tm.VectorUint data type
    fingerprints = [tmap.VectorUint(enc.encode_mol(s)) for s in mols]

    # Initialize the LSH Forest
    lf = tmap.LSHForest(perm)

    # Add the Fingerprints to the LSH Forest and index
    lf.batch_add(fingerprints)
    lf.index()

    nns_per_mol = lf.batch_query(
        [tmap.VectorUint(enc.encode_mol(mols[0]))], 10)

    # Get the coordinates
    x, y, s, t, graph_props = tmap.layout_from_lsh_forest(lf)
    x, y, s, t = np.array(x), np.array(y), np.array(s), np.array(t)

    for nns in nns_per_mol:
        x_nns = x[nns]
        y_nns = y[nns]

    i = 0
    res = []
    xy = list(zip(x, y))
    for x, mol in enumerate(all_mols):
        if not mol:
            res.append(None)
        else:
            res.append(xy[i])
            i += 1

    return res


def tmap_projection(data):
    n = len(data)
    edge_list = []

    max = None
    min = None
    for i, j in zip(range(n), range(n)):
        if i == j:
            continue
        distance = np.linalg.norm(data[i] - data[j])
        edge_list.append([i, j, distance])
        if max is None or distance > max:
            max = distance
        if min is None or distance < min:
            min = distance
    # Normalize
    for edge in edge_list:
        edge[2] = (edge[2] - min) / (max - min)

    x, y, s, t, _ = tmap.layout_from_edge_list(
        n, edge_list, create_mst=False
    )

    return np.array([xy for xy in zip(x, y)])
