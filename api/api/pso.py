from functools import partial
import numpy as np
from flask.views import MethodView
from rdkit import Chem
from mso.objectives.scoring import ScoringFunction
from mso.objectives.mol_functions import qed_score, heavy_atom_count, substructure_match_score, penalize_macrocycles
from mso.optimizer import BasePSOptimizer
from ..schema import PSOArgsSchema, ParticleSchema
from ..utils import cached, catch_time, mol
from ..constants import inference_model, blp, logger
from ..projection import compute_all_projections
from ..registry import models_by_name


def get_scoring_function_from_dict(dictionary):
    name = dictionary['name']
    desirability = dictionary.get('desirability', None)
    weight = dictionary.get('weight', 100)
    kwargs = dictionary.get('additional_args', {})
    func, description, is_mol_func = models_by_name[name]
    if kwargs:
        if name == "distance score":
            target = inference_model.seq_to_emb(kwargs["query"])
            func = partial(func, target=target)
        elif (name == "substructure match") | (name == "substructure exclusion"):
            query = Chem.MolFromSmiles(kwargs["query"])
            func = partial(func, query=query)
        else:
            func = partial(func, **kwargs)
    return ScoringFunction(
        func=func,
        name=name,
        description=description,
        weight=weight,
        desirability=desirability,
        is_mol_func=is_mol_func) if func else None


@blp.route('/pso/')
class PSOAPI(MethodView):

    @blp.arguments(PSOArgsSchema, location='json')
    @blp.response(200, ParticleSchema(many=True))
    @cached
    def post(self, args):
        structure = args.get('structure')
        num_part = args.get('num_part')
        num_swarms = args.get('num_swarms')
        iterations = args.get('iterations')
        v_min = args.get('v_min')
        v_max = args.get('v_max')
        inertia_weight = args.get('inertia_weight')
        phi1 = args.get('phi1')
        phi2 = args.get('phi2')
        phi3 = args.get('phi3')
        objectives = args.get('objectives')

        scoring_functions = [get_scoring_function_from_dict(
            dictionary) for dictionary in objectives]

        logger.info(
            f"Starting PSO from {structure} with {num_part} particles in {num_swarms} swarms for {iterations} iterations")

        hac_desirability = [{"x": 0, "y": 0}, {"x": 5, "y": 0.1}, {"x": 15, "y": 0.9}, {"x": 20, "y": 1.0}, {
            "x": 25, "y": 1.0}, {"x": 30, "y": 0.9, }, {"x": 40, "y": 0.1}, {"x": 45, "y": 0.0}]
        substructure_match = partial(substructure_match_score, query=Chem.MolFromSmiles(
            "CC(=O)N"))  # use partial to define the additional argument (the substructure)
        # invert the resulting score to penalize for a match.
        miss_match_desirability = [{"x": 0, "y": 0}, {"x": 1, "y": 1}]
        # scoring_functions = [
        #     ScoringFunction(heavy_atom_count, "hac",
        #                     desirability=hac_desirability, is_mol_func=True),
        #     ScoringFunction(qed_score, "qed", is_mol_func=True),
        #     ScoringFunction(substructure_match, "miss_match",
        #                     desirability=miss_match_desirability, is_mol_func=True),
        #     ScoringFunction(penalize_macrocycles, "macro", is_mol_func=True)
        # ]

        opt = BasePSOptimizer.from_query(
            init_smiles=structure,
            num_part=num_part,
            num_swarms=num_swarms,
            inference_model=inference_model,
            scoring_functions=scoring_functions,
            v_min=v_min,
            v_max=v_max,
            # inertia_weight=inertia_weight,
            phi1=phi1,
            phi2=phi2,
            phi3=phi3
        )
        # inertia_weight is missing from keyword arguments: https://github.com/jrwnter/mso/blob/master/mso/optimizer.py#L130
        for swarm in opt.swarms:
            swarm.inertia_weight = inertia_weight

        x_stacked = None
        v_stacked = None
        smiles = []
        iteration = []
        instance = []
        swarm_ids = []
        fitness = []
        additional_scores = {}
        for i in range(iterations):
            with catch_time(f'Iteration {i+1}'):
                opt.run(1)
                for swarm_id, swarm in enumerate(opt.swarms):
                    x_stacked = np.concatenate(
                        list(filter(lambda x: x is not None, [x_stacked, swarm.x])), axis=0)
                    v_stacked = np.concatenate(
                        list(filter(lambda x: x is not None, [v_stacked, swarm.v])), axis=0)
                    smiles += swarm.smiles
                    iteration += [i] * swarm.num_part
                    instance += [inst for inst in list(
                        range(swarm_id * num_part, swarm_id * num_part + swarm.num_part))]
                    swarm_ids += [swarm_id] * swarm.num_part
                    fitness += swarm.fitness.tolist()
                    for key in swarm.unscaled_scores.keys():
                        additional_scores[key] = additional_scores.get(
                            key, []) + swarm.unscaled_scores[key].tolist()
                        additional_scores[f'{key}_desirable'] = additional_scores.get(
                            f'{key}_desirable', []) + swarm.desirability_scores[key].tolist()

        # Convert, otherwise "Array must be symmetric" is thrown by MDS
        x_stacked = x_stacked.astype(np.float64).tolist()
        with catch_time('Computing projections'):
            projections, projection = compute_all_projections(
                {'smiles': smiles, 'embedding': x_stacked}, None)

        return [{
            'structure': smiles[i],
            # TODO: We probably don't want to return the embedding...
            # 'embedding': x_stacked[i],
            'projection': {t: projection[t][0][i] for t in projections},
            'properties': {
                **mol.compute_properties(smiles[i]),
                'iteration': iteration[i],
                'instance_num': instance[i],
                'instance': str(instance[i]),
                'swarm_num': swarm_ids[i],
                'swarm': str(swarm_ids[i]),
                # 'v': v_stacked[i].tolist(),
                'fitness': fitness[i],
                **{key: value[i] for key, value in additional_scores.items()},
            }
        } for i, _ in enumerate(smiles)]
