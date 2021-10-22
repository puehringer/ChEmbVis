from ..schema import ProjectionSchema, ProjectionWithModelsArgsSchema
from ..utils import cached, logger
from ..constants import blp
from ..projection import compute_all_projections


@blp.route('/projection/models', methods=["POST"])
@blp.arguments(ProjectionWithModelsArgsSchema, location='json')
@blp.response(200, ProjectionSchema())
@cached
def compute_projection_with_models(args):
    particles = args.get('particles')
    models = args.get('models')

    # Fetch the possible embeddings from the first particle
    possible_embeddings = (particles[0] or {}).get('embedding', {}).keys()

    projections, projection = compute_all_projections(
        None,
        {
            'smiles': [p['structure'] for p in particles],
            **{key: [p.get('embedding', {}).get(key) for p in particles] for key in possible_embeddings}
        },
        {
            'models': models,
            'all_projections': models.keys()
        })

    return {
        'projections': projections,
        'projection': {t: projection[t][1] for t in projections}
    }
