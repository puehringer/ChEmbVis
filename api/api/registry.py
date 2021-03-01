from flask.views import MethodView
from flask import jsonify
from ..utils import cached
from ..constants import blp
from ..registry import model_description


@blp.route('/registry/')
class RegistryAPI(MethodView):

    # TODO: Models
    # @blp.arguments(ProjectionArgsSchema, location='json')
    # @blp.response(200, ProjectionSchema())
    def get(self):
        return jsonify({'objectives': model_description})
        
