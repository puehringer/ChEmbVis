from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.ML.Descriptors import MoleculeDescriptors


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
        **dict(zip(calculator.GetDescriptorNames(), calculator.CalcDescriptors(mol))),
        # 'morgan': list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
    }


def _get_svg_drawer():
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    _options = drawer.drawOptions()
    # Make background transparent
    _options.clearBackground = False
    return drawer


def _get_svg_from_drawer(drawer):
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace("<?xml version='1.0' encoding='iso-8859-1'?>\n", '')


def _draw_mol_to_svg(mol, substructure=None):
    drawer = _get_svg_drawer()
    drawer.DrawMolecule(mol,
                        highlightAtoms=mol.GetSubstructMatch(substructure) if substructure is not None else None)
    return _get_svg_from_drawer(drawer)


def mol_to_svg(mol, substructure=None):
    mol = _string_to_mol(mol)
    substructure = _string_to_mol(substructure)
    if not mol:
        return None
    return _draw_mol_to_svg(mol, substructure=substructure)


def mols_to_mcs_svg(mols):
    mols = list(filter(None, [_string_to_mol(mol) for mol in mols]))
    # TODO: Mols might be null
    res = rdFMCS.FindMCS(mols, matchValences=True, ringMatchesRingOnly=True)
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


def to_fingerprint(fingerprint):
    if fingerprint == 'morgan':
        def fp_getter(m): return AllChem.GetMorganFingerprintAsBitVect(m, 2)
    elif fingerprint == 'daylight':
        def fp_getter(m): return Chem.RDKFingerprint(m)
    else:
        raise ValueError(f'Fingerprint {fingerprint} unknown')
    return fp_getter


def mols_to_fingerprints(mols, fingerprint):
    fp_getter = to_fingerprint(fingerprint)
    mols = [_string_to_mol(mol) for mol in mols]
    return [fp_getter(m) for m in mols]
