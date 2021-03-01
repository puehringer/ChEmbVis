from flask.views import MethodView
from flask import Response
from flask_smorest import abort
from rdkit import DataStructs, Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from ..utils import cached, mol
from ..constants import blp, logger
from ..schema import MoleculesImageArgsSchema, MoleculeImageArgsSchema, MoleculesSubstructureArgsSchema, MoleculesSubstructureSchema, MoleculesTanimotoSchema, MoleculesTanimotoArgsSchema


@blp.route('/image/')
class MoleculeImageAPI(MethodView):

    @blp.arguments(MoleculeImageArgsSchema, location='query')
    @cached
    def get(self, args):
        svg = mol.mol_to_svg(args.get('structure'), args.get('substructure'))
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

        reference_mol = Chem.MolFromSmiles(reference)

        if not reference_mol:
            raise ValueError(
                f'Reference could not be parsed as SMILES: {reference}')

        # Define fingerprint getter
        fp_getter = mol.to_fingerprint(fingerprint)
        reference_fp = fp_getter(reference_mol)

        # TODO: Use DataStructs.BulkTanimotoSimilarity instead
        tanimoto = {}
        for smiles in set(structures):
            m = mol._string_to_mol(smiles)
            if m:
                current_fp = fp_getter(m)
                tanimoto[smiles] = DataStructs.FingerprintSimilarity(
                    reference_fp, current_fp)

        return {'tanimoto': tanimoto}
