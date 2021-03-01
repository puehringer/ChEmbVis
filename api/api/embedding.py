import itertools
from flask import jsonify
from flask_smorest import abort
from flask.views import MethodView
from ..schema import EmbeddingArgsSchema, ParticleSchema
from ..utils import cached, mol, chunkify, catch_time, parallelized
from ..constants import inference_model, blp, logger
from ..projection import compute_all_projections, chembl_umap


@blp.route('/embedding/')
class EmbeddingAPI(MethodView):

    @cached
    def get(self):
        # TODO: Move to different API
        if not chembl_umap.projected_umap or not chembl_umap.projected_pca:
            abort(404)
        structures = list(map(lambda p: p.get('structure'), chembl_umap.projected_umap))
        umap_projections = list(map(lambda p: p.get('projection'), chembl_umap.projected_umap))
        pca_projections = list(map(lambda p: p.get('projection'), chembl_umap.projected_pca))
        properties = list(parallelized(mol.compute_properties, structures))

        return jsonify([{
            'structure': structure,
            'projection': {
                'chembl_umap': umap_projections[i],
                'chembl_pca': pca_projections[i]
            },
            'properties': properties[i]
        } for i, structure in enumerate(structures)])

    @blp.arguments(EmbeddingArgsSchema, location='json')
    @blp.response(200, ParticleSchema(many=True))
    @cached
    def post(self, args):
        structures = args.get('structures')

        with catch_time(f'Filtering {len(structures)} structures'):
            structures = [s for s in structures if mol.is_valid_mol(s)]

        # Add "chunkify" to avoid OOM error
        with catch_time(f'Computing embedding for {len(structures)} structures'):
            # smiles_embedding = inference_model.seq_to_emb(structures[:1000]).tolist()
            smiles_embedding = list(itertools.chain.from_iterable(inference_model.seq_to_emb(chunk).tolist() for chunk in chunkify(structures, 512)))
            # TODO: This is actually not required...
            # decoded_smiles_list = list(itertools.chain.from_iterable([inference_model.emb_to_seq(chunk) for chunk in smiles_embedding_chunks]))

        with catch_time(f'Computing projections for {len(structures)} structures'):
            projections, projection = compute_all_projections(
                {'smiles': structures, 'embedding': smiles_embedding}, None)

        return [{
            'structure': structure,
            # 'decoded': decoded_smiles_list[i],
            # 'embedding': smiles_embedding[i],
            'projection': {t: projection[t][0][i] for t in projections},
            'properties': mol.compute_properties(structure)
        } for i, structure in enumerate(structures)]
