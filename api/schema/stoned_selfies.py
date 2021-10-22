import marshmallow as ma
from . import ParticleSchema

class StonedSelfiesArgsSchema(ma.Schema):
    structure = ma.fields.String(required=True)
    substructure = ma.fields.String(required=False)
    random_samples = ma.fields.Integer(required=False, missing=1000)
    max_mutations = ma.fields.Integer(required=False, missing=5)


# class NeighborhoodSamplingArgsSchema(ma.Schema):
#     structure = ma.fields.String(required=True)
#     samples = ma.fields.Integer(strict=True, missing=10, validate=ma.validate.Range(min=2))
#     method = ma.fields.String(required=False, missing='chembl_pca', validate=ma.validate.OneOf(['chembl_pca', 'random_orthogonal']))
#     scale = ma.fields.Float(missing=1, validate=ma.validate.Range(min=0))


class StonedSelfiesSchema(ParticleSchema):
    structures = ma.fields.List(ma.fields.String())