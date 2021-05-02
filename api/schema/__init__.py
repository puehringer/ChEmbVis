import marshmallow as ma


class EmbeddingArgsSchema(ma.Schema):
    structures = ma.fields.List(
        ma.fields.String(), required=True, validate=ma.validate.Length(min=1))
    include_embedding = ma.fields.Boolean(missing=False)


class ParticleSchema(ma.Schema):
    structure = ma.fields.String(required=True)
    # decoded = ma.fields.String(required=False)
    embedding = ma.fields.List(ma.fields.Float(), required=False)
    projection = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(
        ma.fields.Float()), required=False)
    properties = ma.fields.Raw(required=False)


from .projection import ProjectionSchema, ProjectionArgsSchema
from .pso import PSOArgsSchema
from .interpolation import InterpolationArgsSchema, InterpolatedParticleSchema, NeighborhoodSamplingArgsSchema
from .molecule import *
from .mmp import *

