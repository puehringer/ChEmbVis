from flask.views import MethodView
from flask import Response, jsonify
from flask_smorest import abort
import requests
from rdkit import DataStructs, Chem
from rdkit.Chem import TemplateAlign
from rdkit.Chem.Scaffolds import MurckoScaffold
from ..utils import cached, mol, parallelized
from ..constants import blp
from ..models import embedding_models
from ..schema import MoleculesImageArgsSchema, MoleculeImageArgsSchema, MoleculesSubstructureArgsSchema, MoleculesSubstructureSchema, MoleculesTanimotoSchema, MoleculesTanimotoArgsSchema


@blp.route('/image/')
class MoleculeImageAPI(MethodView):

    @blp.arguments(MoleculeImageArgsSchema, location='query')
    @cached
    def get(self, args):
        substructure = args.get('substructure')
        align_mol = mol._string_to_mol(args.get('align'))
        structure_mol = mol._string_to_mol(args.get('structure'))
        if align_mol and structure_mol:
            TemplateAlign.rdDepictor.Compute2DCoords(structure_mol)
            # Find the maximum common substructure for aligning purposes
            mcs = mol.mols_to_mcs([structure_mol, align_mol]).queryMol
            if mcs:
                TemplateAlign.rdDepictor.Compute2DCoords(mcs)
                TemplateAlign.AlignMolToTemplate2D(structure_mol, mcs, clearConfs=True)
                # Enable this to show substructore highlights of alignment
                # substructure = mcs
        svg = mol.mol_to_svg(structure_mol, substructure)
        if not svg:
            abort(404)
        return Response(svg, mimetype='image/svg+xml')

    @blp.arguments(MoleculesImageArgsSchema, location='json')
    @cached
    def post(self, args):
        structures = args.get('structures')
        method = args.get('method')
        substructure = args.get('substructure')

        mols = list(filter(None, [mol._string_to_mol(s)
                                  for s in set(structures)]))

        if not mols:
            abort(404)

        if method == 'auto':
            method = 'mcs' if len(mols) > 2 else (
                'similarity' if len(mols) == 2 else 'single')

        if method == 'single':
            svg = mol.mol_to_svg(mols[0], substructure)
        elif method == 'murcko':
            # https://www.rdkit.org/docs/GettingStartedInPython.html#murcko-decomposition
            # MurckoScaffold.MakeScaffoldGeneric(core)
            svg = mol.mol_to_svg(MurckoScaffold.GetScaffoldForMol(mols[0]))
        elif method == 'mcs':
            # https://www.rdkit.org/docs/GettingStartedInPython.html#maximum-common-substructure
            svg = mol.mols_to_mcs_svg(mols)
        elif method == 'similarity':
            svg = mol.mols_to_similarity_svg(mols)

        if not svg:
            abort(404)
        return Response(svg, mimetype='image/svg+xml')


@blp.route('/mol/substructures/')
class MoleculeAPI(MethodView):

    @blp.arguments(MoleculesSubstructureArgsSchema, location='json')
    @blp.response(200, MoleculesSubstructureSchema())
    @cached
    def post(self, args):
        structures = args.get('structures')
        smarts = args.get('smarts')

        smarts_substructure = Chem.MolFromSmarts(smarts)
        smiles_substructure = Chem.MolFromSmiles(smarts)

        if not smarts_substructure and not smiles_substructure:
            raise ValueError(
                f'Input could not be parsed as SMARTS or SMILES: {smarts}')

        validity = {}
        counts = {}
        for smiles in set(structures):
            m = mol._string_to_mol(smiles)
            validity[smiles] = bool(m and (smarts_substructure and m.HasSubstructMatch(
                smarts_substructure) or (smiles_substructure and m.HasSubstructMatch(smiles_substructure))))

            counts[smiles] = max(0, bool(m and smarts_substructure) and len(m.GetSubstructMatches(
                smarts_substructure)), bool(m and smiles_substructure) and len(m.GetSubstructMatches(smiles_substructure)))

        return {'validity': validity, 'counts': counts}


@blp.route('/mol/tanimoto/')
class TanimotoAPI(MethodView):

    @blp.arguments(MoleculesTanimotoArgsSchema, location='json')
    @blp.response(200, MoleculesTanimotoSchema())
    @cached
    def post(self, args):
        structures = args.get('structures')
        reference = args.get('reference')
        fingerprint = args.get('fingerprint')

        if fingerprint == 'cddd':
            return jsonify({'tanimoto': requests.post('http://api_umap:5000/api/similarity/', json={'reference': embedding_models['cddd'].encode_single(reference) }).json()})

        reference_mol = Chem.MolFromSmiles(reference)

        if not reference_mol:
            raise ValueError(
                f'Reference could not be parsed as SMILES: {reference}')

        # Define fingerprint getter
        fp_getter = mol.to_fingerprint(fingerprint)
        reference_fp = fp_getter(reference_mol)

        # Generate mol and fingerprint
        def _get_fp(smiles: str):
            m = mol._string_to_mol(smiles)
            if not m:
                return None
            return (smiles, fp_getter(m))
        # Compute fingerprints
        smiles_to_fingerprints = dict(filter(None, parallelized(_get_fp, set(structures))))
        # Compute similarities in bulk
        similarities = DataStructs.BulkTanimotoSimilarity(reference_fp, list(smiles_to_fingerprints.values()))
        # Return smiles -> similarities object
        return {'tanimoto': {s: similarities[i] for i, s in enumerate(smiles_to_fingerprints.keys())}}
