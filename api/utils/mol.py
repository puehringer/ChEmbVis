from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, Mol
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.AtomPairs.Sheridan import GetBPFingerprint, GetBTFingerprint
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem import MolToSmiles
import numpy as np
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_valid_mol(mol):
    return _string_to_mol(mol) is not None


def _string_to_mol(mol):
    if isinstance(mol, Chem.rdchem.Mol):
        return mol
    elif isinstance(mol, str):
        return Chem.MolFromSmiles(mol)
    return None


def compute_properties(mol):
    mol = _string_to_mol(mol)
    if not mol:
        return {'valid': False}

    # all_available_descriptors = [x[0] for x in Chem.Descriptors._descList]
    all_available_descriptors = ['MolLogP', 'MolWt', 'MolMR', 'TPSA', 'qed']
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(
        all_available_descriptors)

    return {
        'valid': True,
        'scaffold': _generate_scaffold(MolToSmiles(mol)),
        **dict(zip(calculator.GetDescriptorNames(), calculator.CalcDescriptors(mol)))
    }


def _generate_scaffold(smiles: str, include_chirality: bool = False) -> str:
    """Compute the Bemis-Murcko scaffold for a SMILES string.
    Bemis-Murcko scaffolds are described in DOI: 10.1021/jm9602928.
    They are essentially that part of the molecule consisting of
    rings and the linker atoms between them.
    Paramters
    ---------
    smiles: str
        SMILES
    include_chirality: bool, default False
        Whether to include chirality in scaffolds or not.
    Returns
    -------
    str
        The MurckScaffold SMILES from the original SMILES
    References
    ----------
    .. [1] Bemis, Guy W., and Mark A. Murcko. "The properties of known drugs.
        1. Molecular frameworks." Journal of medicinal chemistry 39.15 (1996): 2887-2893.
    Note
    ----
    This function requires RDKit to be installed.
    """
    from rdkit import Chem
    from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles

    mol = Chem.MolFromSmiles(smiles)
    scaffold = MurckoScaffoldSmiles(mol=mol, includeChirality=include_chirality)
    return scaffold


def _get_svg_drawer():
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    _options = drawer.drawOptions()
    # Make background transparent
    _options.clearBackground = False
    return drawer


def _get_svg_from_drawer(drawer):
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace("<?xml version='1.0' encoding='iso-8859-1'?>\n", '')


def _draw_mol_to_svg(mol, substructure=None, highlight_atoms=None, highlight_bonds=None, highlight_atom_colors=None, highlight_bond_colors=None):
    drawer = _get_svg_drawer()
    drawer.DrawMolecule(mol,
                        highlightAtoms=highlight_atoms if highlight_atoms else (mol.GetSubstructMatch(substructure) if substructure is not None else None),
                        highlightBonds=highlight_bonds,
                        highlightAtomColors=highlight_atom_colors,
                        highlightBondColors=highlight_bond_colors)
    return _get_svg_from_drawer(drawer)


def mol_to_svg(mol, substructure=None, highlight_atoms=None, highlight_bonds=None, highlight_atom_colors=None, highlight_bond_colors=None):
    mol = _string_to_mol(mol)
    substructure = _string_to_mol(substructure)
    if not mol:
        return None
    return _draw_mol_to_svg(mol, substructure=substructure, highlight_atoms=highlight_atoms, highlight_bonds=highlight_bonds, highlight_atom_colors=highlight_atom_colors, highlight_bond_colors=highlight_bond_colors)


def mols_to_mcs(mols):
    mols = list(filter(None, [_string_to_mol(mol) for mol in mols]))
    return rdFMCS.FindMCS(mols, matchValences=True, ringMatchesRingOnly=True, completeRingsOnly=True)


def mols_to_mcs_svg(mols):
    res = mols_to_mcs(mols)
    return _draw_mol_to_svg(res.queryMol)


def mols_to_similarity_svg(mols):
    mols = [_string_to_mol(mol) for mol in mols]
    drawer = _get_svg_drawer()
    ref, probe = mols[0], mols[1] if len(mols) > 1 else mols[0]
    # Compute similarity between two molecules: https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/Draw/SimilarityMaps.py
    _, maxWeight = SimilarityMaps.GetSimilarityMapForFingerprint(ref, probe,
                                                                 lambda m, i: SimilarityMaps.GetMorganFingerprint(
                                                                     m, i, radius=2, fpType='bv'),
                                                                 draw2d=drawer)
    return _get_svg_from_drawer(drawer)


class _FingerprintCalculator:
    ''' Calculate the fingerprint for a molecule, given the fingerprint type
    Parameters: 
        fingerprint (string)            :Fingerprint type  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6)  
    Returns:
        function accepting a MOL and returning the corresponding fingerprint.
    '''

    def get_fingerprint(self, fingerprint: str):
        method_name = 'get_' + fingerprint.upper()
        method = getattr(self, method_name)
        if method is None:
            raise Exception(f'{fingerprint} is not a supported fingerprint type.')
        return method

    def get_AP(self, mol: Mol):
        return AllChem.GetAtomPairFingerprint(mol, maxLength=10)

    def get_PHCO(self, mol: Mol):
        return Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)

    def get_BPF(self, mol: Mol):
        return GetBPFingerprint(mol)

    def get_BTF(self, mol: Mol):
        return GetBTFingerprint(mol)

    def get_PATH(self, mol: Mol):
        return AllChem.RDKFingerprint(mol)

    def get_ECFP2(self, mol: Mol):
        return AllChem.GetMorganFingerprintAsBitVect(mol, 1)

    def get_ECFP4(self, mol: Mol):
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2)

    def get_ECFP6(self, mol: Mol):
        return AllChem.GetMorganFingerprintAsBitVect(mol, 3)

    def get_FCFP4(self, mol: Mol):
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True)

    def get_FCFP6(self, mol: Mol):
        return AllChem.GetMorganFingerprintAsBitVect(mol, 3, useFeatures=True)


def to_fingerprint(fingerprint: str):
    ''' Fingerprint getter method. Fingerprint is returned after using object of 
        class '_FingerprintCalculator'
        
    Parameters: 
        fingerprint (string)            :Fingerprint type  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6)  
    Returns:
        function accepting a MOL and returning the corresponding fingerprint.
        
    '''
    return _FingerprintCalculator().get_fingerprint(fingerprint)


def mols_to_fingerprints(mols, fingerprint: str):
    fp_getter = to_fingerprint(fingerprint)
    mols = [_string_to_mol(mol) or _string_to_mol('*') for mol in mols]
    return [fp_getter(m) for m in mols]

# Preprocessing pipeline from https://github.com/jrwnter/cddd/blob/master/cddd/preprocessing.py
REMOVER = SaltRemover()
ORGANIC_ATOM_SET = set([5, 6, 7, 8, 9, 15, 16, 17, 35, 53])

def randomize_smile(sml):
    """Function that randomizes a SMILES sequnce. This was adapted from the
    implemetation of E. Bjerrum 2017, SMILES Enumeration as Data Augmentation
    for Neural Network Modeling of Molecules.
    Args:
        sml: SMILES sequnce to randomize.
    Return:
        randomized SMILES sequnce or
        nan if SMILES is not interpretable.
    """
    try:
        m = Chem.MolFromSmiles(sml)
        ans = list(range(m.GetNumAtoms()))
        np.random.shuffle(ans)
        nm = Chem.RenumberAtoms(m, ans)
        return Chem.MolToSmiles(nm, canonical=False)
    except:
        return float('nan')

def keep_largest_fragment(sml):
    """Function that returns the SMILES sequence of the largest fragment for a input
    SMILES sequnce.
    Args:
        sml: SMILES sequence.
    Returns:
        canonical SMILES sequnce of the largest fragment.
    """
    mol_frags = Chem.GetMolFrags(Chem.MolFromSmiles(sml), asMols=True)
    largest_mol = None
    largest_mol_size = 0
    for mol in mol_frags:
        size = mol.GetNumAtoms()
        if size > largest_mol_size:
            largest_mol = mol
            largest_mol_size = size
    return Chem.MolToSmiles(largest_mol)

def remove_salt_stereo(sml, remover):
    """Function that strips salts and removes stereochemistry information from a SMILES.
    Args:
        sml: SMILES sequence.
        remover: RDKit's SaltRemover object.
    Returns:
        canonical SMILES sequnce without salts and stereochemistry information.
    """
    try:
        sml = Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(sml),
                                                dontRemoveEverything=True),
                               isomericSmiles=False)
        if "." in sml:
            sml = keep_largest_fragment(sml)
    except:
        sml = np.float("nan")
    return(sml)

def filter_smiles(sml):
    try:
        m = Chem.MolFromSmiles(sml)
        logp = Descriptors.MolLogP(m)
        mol_weight = Descriptors.MolWt(m)
        num_heavy_atoms = Descriptors.HeavyAtomCount(m)
        atom_num_list = [atom.GetAtomicNum() for atom in m.GetAtoms()]
        is_organic = set(atom_num_list) <= ORGANIC_ATOM_SET
        if ((logp > -5) & (logp < 7) &
            (mol_weight > 12) & (mol_weight < 600) &
            (num_heavy_atoms > 3) & (num_heavy_atoms < 50) &
            is_organic ):
            return Chem.MolToSmiles(m)
        else:
            return float('nan')
    except:
        return float('nan')

def preprocess_smiles(sml):
    """Function that preprocesses a SMILES string such that it is in the same format as
    the translation model was trained on. It removes salts and stereochemistry from the
    SMILES sequnce. If the sequnce correspond to an inorganic molecule or cannot be
    interpreted by RDKit nan is returned.
    Args:
        sml: SMILES sequence.
    Returns:
        preprocessd SMILES sequnces or nan.
    """
    new_sml = remove_salt_stereo(sml, REMOVER)
    new_sml = filter_smiles(new_sml)
    return new_sml