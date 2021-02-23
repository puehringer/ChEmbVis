from flask.views import MethodView
import numpy as np
from ..schema import InterpolationArgsSchema, ParticleSchema
from ..utils import cached, mol
from ..constants import inference_model, blp
from ..projection import compute_all_projections


@blp.route('/interpolation/')
class InterpolationAPI(MethodView):

    @blp.arguments(InterpolationArgsSchema, location='json')
    @blp.response(200, ParticleSchema(many=True))
    @cached
    def post(self, args):
        max_samples = args.get('maxSamples')
        structures = args.get('structures')
        # Compute embedding for all structures
        embedded_structures = inference_model.seq_to_emb(structures)
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
        interpolated = inference_model.emb_to_seq(
            np.array(embedded_interpolated))

        # Store already seen structures
        seen = set()
        structures = []
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
                'embedding': embedded,
                'scaffold': scaffold,
                'properties': mol.compute_properties(structure)
            })

        projections, projection = compute_all_projections(
            {'smiles': list(map(lambda s: s['structure'], structures)), 'embedding': list(map(lambda s: s['embedding'], structures))}, None)

        return [{**o, 'projection': {t: projection[t][0][i] for t in projections}} for i, o in enumerate(structures)]
