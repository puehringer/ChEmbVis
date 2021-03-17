from flask.views import MethodView
from subprocess import PIPE, Popen
import io
from ..schema import MMPArgsSchema, MMPSchema
from ..utils import cached
from ..constants import blp, logger
from ..projection import compute_all_projections


@blp.route('/mmp/')
class MMPAPI(MethodView):

    @blp.arguments(MMPArgsSchema, location='json')
    @blp.response(200, MMPSchema())
    @cached
    def post(self, args):
        structure = args.get('structure')
        min_variable_size = args.get('min_variable_size')
        max_variable_size = args.get('max_variable_size')
        min_constant_size = args.get('min_constant_size')
        min_radius = args.get('min_radius')
        min_pairs = args.get('min_pairs')
        substructure = args.get('substructure')

        command_str = f"mmpdb transform --smiles '{structure}' --min-variable-size {min_variable_size} --max-variable-size {max_variable_size} --min-constant-size {min_constant_size} --min-radius {min_radius} --min-pairs {min_pairs} --substructure '{substructure}' /_shared/chembl.mmpdb"
        logger.info(f'Invoking {command_str}')

        command = Popen(command_str, stdout=PIPE, stderr=PIPE, shell=True)
        errors = list([line.decode() for line in command.stderr])
        output = list([line.decode() for line in command.stdout])
        if not output and errors:
            raise ValueError(f'mmpdb returned errors: {", ".join(errors)}')

        # Hack: parse the second argument (the smiles) for each line
        structures = list(set([line.strip().split('\t')[1]
                               for line in output[1:]]))

        return {
            'structures': structures
        }
