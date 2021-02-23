from flask import jsonify
from flask_smorest import abort
from flask.views import MethodView
from ..schema import EmbeddingArgsSchema, ParticleSchema
from ..utils import cached, mol
from ..constants import inference_model, blp
from ..projection import compute_all_projections, chembl_umap


@blp.route('/embedding/')
class EmbeddingAPI(MethodView):

    @cached
    def get(self):
        # TODO: Move to different API
        if not chembl_umap.projected_umap:
            abort(404)
        # return jsonify(chembl_umap.projected_umap)
        structures = list(map(lambda p: p.get('structure'),
                              chembl_umap.projected_umap))
        projections = list(
            map(lambda p: p.get('projection'), chembl_umap.projected_umap))
        return jsonify([{
            'structure': structure,
            'projection': projections[i],
            'properties': mol.compute_properties(structure)
        } for i, structure in enumerate(structures[:15_000]) if mol.is_valid_mol(structure)])

    @blp.arguments(EmbeddingArgsSchema, location='json')
    @blp.response(200, ParticleSchema(many=True))
    @cached
    def post(self, args):
        structures = args.get('structures')
        additional = args.get('additional', {})
        smiles_embedding = inference_model.seq_to_emb(
            args.get('structures'))
        decoded_smiles_list = inference_model.emb_to_seq(smiles_embedding)

        smiles_embedding = smiles_embedding.tolist()

        projections, projection = compute_all_projections(
            {'smiles': structures, 'embedding': smiles_embedding}, None)

        return [{
            'structure': structure,
            'decoded': decoded_smiles_list[i],
            'embedding': smiles_embedding[i],
            'projection': {t: projection[t][0][i] for t in projections},
            'properties': {
                **mol.compute_properties(structure),
                **{key: value[i] for key, value in additional.items()},
            }
        } for i, structure in enumerate(structures) if mol.is_valid_mol(structure)]
