from flask.views import MethodView
from ..schema import ProjectionSchema, ProjectionArgsSchema
from ..utils import cached
from ..constants import blp
from ..projection import compute_all_projections


@blp.route('/projection/')
class ProjectionAPI(MethodView):

    @blp.arguments(ProjectionArgsSchema, location='json')
    @blp.response(200, ProjectionSchema())
    @cached
    def post(self, args):
        structures = args.get('structures')
        embedding = args.get('embedding')
        additional_structures = args.get('additional_structures')
        additional_embedding = args.get('additional_embedding')
        projections, projection = compute_all_projections({'smiles': structures, 'embedding': embedding}, {
                                                          'smiles': additional_structures, 'embedding': additional_embedding})

        return {
            'projections': projections,
            'projection': {t: projection[t][0] for t in projections},
            'additional': {t: projection[t][1] if projection[t][1] is not None else None for t in projections},
        }
