import numpy as np
from scipy.stats import ortho_group
from rdkit.Chem import TemplateAlign
from ..schema import InterpolationArgsSchema, NeighborhoodSamplingArgsSchema, ServerCollectionSchema
from ..utils import cached, mol, logger
from ..constants import blp
from ..models import embedding_models
from ..projection import compute_all_projections
from ..projection.chembl_umap import get_projected_pca_components


@blp.route('/sampling/', methods=['POST'])
@blp.arguments(NeighborhoodSamplingArgsSchema, location='json')
@blp.response(200, ServerCollectionSchema)
@cached
def neighborhood_sampling(args):
    samples = args.get('samples')
    structure = args.get('structure')
    method = args.get('method')
    scale = args.get('scale')

    logger.info(
        f'Starting {method} neighborhood sampling with {samples} samples at scale {scale} from structure {structure}')

    # Compute embedding for structure
    embedded_structure = embedding_models['cddd'].encode_single(structure)

    dim = len(embedded_structure)
    # TODO: Find proper vectors and factor
    if method == 'chembl_pca':
        vectors = get_projected_pca_components()
    elif method == 'random_orthogonal':
        vectors = ortho_group.rvs(dim)[:2]
    # vectors = np.array([np.ones((512,)), np.zeros((512,))]) / samples
    # vectors = (2 * np.random.rand(2, 512) - 1) / samples

    vectors *= scale

    # Check length of vector (should be 1 * factor)
    # logger.info(np.linalg.norm(vectors, axis=1))

    # Create grid ranges
    scaler = np.arange(0, samples) - samples // 2
    # Create pairs of grid coordinates
    pairs = np.array(np.meshgrid(scaler, scaler)).T.reshape(-1, 2)
    # Swap X and Y axis as [[-5, -5], [-5, -4], ...] is currently in pairs
    pairs = pairs[:, [1, 0]]
    # Invert the Y axis to get a proper "grid", i.e. [[-5, 5], [-4, 5], ...]
    pairs[:, 1] *= -1
    # Create grid embedding values: 20x2 @ 2x512 = 20x512
    result = (pairs @ vectors) + embedded_structure

    # Create grid. TODO .reshape(samples, samples)
    smiles = embedding_models['cddd'].decode(result)
    # Ensure that center structure is the reference structure.
    # This is required as the decode step is not 100% correct.
    smiles[len(smiles) // 2] = structure

    structures = [{
        'structure': structure,
        'embedding': { 'cddd': result[i].tolist() },
        'properties': mol.compute_properties(structure)
    } for i, structure in enumerate(smiles)]

    projections, projection = compute_all_projections({'smiles': list(map(
        lambda s: s['structure'], structures)), 'cddd': list(map(lambda s: s['embedding']['cddd'], structures))}, None)

    return {
        "data": [{**o, 'projection': {t: projection[t][0][i] for t in projections}} for i, o in enumerate(structures)],
        "projections": projections
    }


@blp.route('/interpolation/', methods=["POST"])
@blp.arguments(InterpolationArgsSchema, location='json')
@blp.response(200, ServerCollectionSchema)
@cached
def interpolation(args):
    max_samples = args.get('maxSamples')
    structures = args.get('structures')
    # Compute embedding for all structures
    embedded_structures = embedding_models['cddd'].encode(structures)
    embedded_interpolated = []
    real_structure = []
    # For all consecutive pairs
    for embedded_origin, embedded_target, origin, target in zip(embedded_structures, embedded_structures[1:], structures, structures[1:]):
        # Compute delta
        delta = (embedded_target - embedded_origin) / (max_samples - 1)
        # Linearly interpolate between origin and target
        embedded_interpolated += [(embedded_origin + i * delta).tolist()
                                  for i in range(max_samples)]
        # Save which structures were scaffolds
        real_structure += ([False] + [False] * (max_samples - 2) + [True])
    # Set the first structure to scaffold
    real_structure[0] = True

    # Convert to smiles
    interpolated = embedding_models['cddd'].decode(np.array(embedded_interpolated))

    # Store already seen structures
    seen = set()
    structures = []
    index = 0
    for structure, next_structure, embedded, scaffold in zip(interpolated, [*interpolated[1:], None], embedded_interpolated, real_structure):
        # If we have one of the scaffold structures, clear the seen structures
        if scaffold:
            seen.clear()
        # If a structure was already seen, or it is no scaffold and equals the next structure, skip it
        if structure in seen or (not scaffold and structure == next_structure):
            continue
        # Skip invalid mols
        if not mol.is_valid_mol(structure):
            continue
        # Otherwise, save the structure
        seen.add(structure)
        structures.append({
            'structure': structure,
            'embedding': { 'cddd': embedded },
            'scaffold': scaffold,
            'properties': {**mol.compute_properties(structure), "index": index}
        })
        index += 1

    # Find the maximum common substructure for aligning purposes
    mcs = mol.mols_to_mcs(
        map(lambda m: m.get('structure'), structures)).queryMol
    if mcs:
        # Compute 2D coords for aligning
        TemplateAlign.rdDepictor.Compute2DCoords(mcs)

    previous_mol = None
    for structure in structures:
        smiles = structure.get('structure')
        current_mol = mol._string_to_mol(smiles)
        # Compute 2D coords for each interpolated structure
        TemplateAlign.rdDepictor.Compute2DCoords(current_mol)
        if mcs:
            # And align it according to the MCS
            TemplateAlign.AlignMolToTemplate2D(
                current_mol, mcs, clearConfs=True)
        # Compute the highlights from the previous structure
        highlight_atoms = None
        highlight_atom_colors = None
        highlight_bonds = None
        highlight_bond_colors = None
        structure['images'] = [None, None]
        if previous_mol:
            # Compute maximum common substructure for current pair
            pair_mcs = mol.mols_to_mcs(
                [current_mol, previous_mol]).queryMol
            if pair_mcs:
                # Substract all atoms which are not part of the substructure match
                highlight_atoms = set(range(
                    len(current_mol.GetAtoms()))) - set(current_mol.GetSubstructMatch(pair_mcs))
                # highlight_atom_colors = {id: (0.8, 1.0, 0.8) for id in highlight_atoms}
                # Get the highlighted bonds, i.e. the bonds leading to or coming from a highlighted atom
                highlight_bonds = [bond.GetIdx() for bond in current_mol.GetBonds(
                ) if bond.GetBeginAtomIdx() in highlight_atoms or bond.GetEndAtomIdx() in highlight_atoms]
                # highlight_bond_colors = {id: (0.8, 0.9, 0.8) for id in highlight_bonds}

            structure['images'][1] = mol.mols_to_similarity_svg(
                [current_mol, previous_mol])
        # Store the image explicitly
        structure['images'][0] = mol.mol_to_svg(current_mol, highlight_atoms=highlight_atoms, highlight_bonds=highlight_bonds,
                                                highlight_atom_colors=highlight_atom_colors, highlight_bond_colors=highlight_bond_colors)
        # Store the previous mol for a difference highlighting
        previous_mol = current_mol

    projections, projection = compute_all_projections(
        {'smiles': list(map(lambda s: s['structure'], structures)), 'cddd': list(map(lambda s: s['embedding']['cddd'], structures))}, None)

    return {
        "data": [{**o, 'projection': {t: projection[t][0][i] for t in projections}} for i, o in enumerate(structures)],
        "projections": projections
    }
