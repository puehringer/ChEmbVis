from abc import ABC, abstractmethod
from typing import List, Dict
from ..utils import chunkify, ensure_np_array, mol
from ..constants import logger, get_cddd_model, get_vae_model, get_molbert_model, get_graph_model
import itertools
import numpy as np


class Embedding(ABC):
    def _method_implemented(self, method: str):
        this_method = getattr(self, method, None)
        base_method = getattr(Embedding, method, None)
        return this_method is not None and this_method.__func__ is not base_method

    @property
    @abstractmethod
    def distance_metric(self):
        pass

    @property
    def can_encode(self):
        return self._method_implemented('encode')

    @property
    def can_decode(self):
        return self._method_implemented('decode')

    @property
    def can_interpolate(self):
        return self._method_implemented('interpolate')

    @property
    def use_raw_smiles(self):
        return False

    def encode(self, smiles: List[str]) -> np.ndarray:
        pass

    def encode_single(self, smiles: str) -> List[float]:
        return self.encode([smiles])[0]

    def decode(self, embedding: List[List[float]]) -> List[str]:
        pass

    def decode_single(self, embedding: List[float]) -> str:
        return self.decode([embedding])[0]


class CDDDEmbedding(Embedding):
    @property
    def model(self):
        return get_cddd_model()

    @property
    def distance_metric(self):
        return 'euclidean'

    def encode(self, smiles: List[str]) -> np.ndarray:
        # Add "chunkify" to avoid OOM error
        return np.array(list(itertools.chain.from_iterable(
            self.model.seq_to_emb(chunk) for chunk in chunkify(smiles, 512))))

    def decode(self, embedding: List[List[float]]) -> List[str]:
        return list(itertools.chain.from_iterable([self.model.emb_to_seq(chunk) for chunk in chunkify(embedding, 512)]))


class GraphEmbedding(Embedding):
    @property
    def distance_metric(self):
        return 'euclidean'

    @property
    def use_raw_smiles(self):
        return True

    def encode(self, smiles: List[str]) -> np.ndarray:
        return get_graph_model()(smiles)


class ChemNetEmbedding(Embedding):
    @property
    def model(self):
        from fcd import load_ref_model
        return load_ref_model()

    @property
    def distance_metric(self):
        return 'euclidean'

    def encode(self, smiles: List[str]) -> np.ndarray:
        from fcd import get_predictions
        return get_predictions(self.model, smiles)


class MOLBertEmbedding(Embedding):
    @property
    def model(self):
        return get_molbert_model()

    @property
    def distance_metric(self):
        return 'euclidean'

    def encode(self, smiles: List[str]) -> np.ndarray:
        return self.model.transform(smiles)[0]


class VAEEmbedding(Embedding):
    @property
    def model(self):
        return get_vae_model()

    @property
    def distance_metric(self):
        return 'cosine'

    def encode(self, smiles: List[str]) -> np.ndarray:
        # Example from github:
        # smiles_1 = mu.canon_smiles('CSCC(=O)NNC(=O)c1c(C)oc(C)c1C')
        # X_1 = vae.smiles_to_hot(smiles_1,canonize_smiles=True)
        # z_1 = vae.encode(X_1)
        # X_r= vae.decode(z_1)

        def transform(smiles):
            from chemvae import mol_utils as mu
            smiles = [mu.canon_smiles(s) for s in smiles]

            filtered_smiles = []
            for smile in smiles:
                converted = mu.smiles_to_hot_filter([smile], self.model.char_indices)
                if len(converted) == 1:
                    filtered_smiles.append(smile)
                else:
                    # TODO: Remove dummy
                    logger.info(f'VAE could not handle a character in {smile}, setting default...')
                    filtered_smiles.append("C1CCCCC1")

            return self.model.encode(self.model.smiles_to_hot(filtered_smiles, canonize_smiles=False, check_smiles=False))

        # import numpy as np
        # vae_embedding = np.zeros((len(structures), 200)).tolist()
        return transform(smiles)

    def decode(self, embedding: List[List[float]]) -> List[str]:
        return self.model.hot_to_smiles(self.model.decode(np.array(embedding)), strip=True)


class FingerprintEmbedding(Embedding):

    def __init__(self, fp: str):
        super().__init__()
        self.fp = fp

    @property
    def distance_metric(self):
        return 'jaccard'

    def encode(self, smiles: List[str]) -> np.ndarray:
        return ensure_np_array(mol.mols_to_fingerprints(smiles, self.fp))


class MAP4Embedding(Embedding):

    def __init__(self):
        super().__init__()
        import tmap as tm
        from map4 import MAP4Calculator

        dim = 1024
        self.MAP4 = MAP4Calculator(dimensions=dim)
        self.ENC = tm.Minhash(dim)

    def distance_metric(self, x, y):
        return self.ENC.get_distance(x, y)

    def encode(self, smiles: List[str]) -> np.ndarray:
        mols = [mol._string_to_mol(m) or mol._string_to_mol('*') for m in mols]
        return ensure_np_array(self.MAP4.calculate_many(mols))


embedding_models: Dict[str, Embedding] = {
    # 'cddd': CDDDEmbedding(),
    # 'graph': GraphEmbedding(),
    # 'chemnet': ChemNetEmbedding(),
    # 'molbert': MOLBertEmbedding(),
    # 'vae': VAEEmbedding(),
    'maccs': FingerprintEmbedding('maccs'),
    # 'ecfp0': FingerprintEmbedding('ecfp0'),
    # 'ecfp2': FingerprintEmbedding('ecfp2'),
    'ecfp4': FingerprintEmbedding('ecfp4'),
    # 'ecfp6': FingerprintEmbedding('ecfp6'),
    # 'rdkit': FingerprintEmbedding('path'),
    # 'mhfp2': FingerprintEmbedding('mhfp2'),
    # 'mhfp4': FingerprintEmbedding('mhfp4'),
    # 'mhfp6': FingerprintEmbedding('mhfp6'),
    # 'ecfp2_256': FingerprintEmbedding('ecfp2_256'),
    # 'ecfp4_256': FingerprintEmbedding('ecfp4_256'),
    # 'ecfp6_256': FingerprintEmbedding('ecfp6_256'),
    # 'map4': MAP4Embedding()
}