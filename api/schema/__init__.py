import marshmallow as ma


class EmbeddingArgsSchema(ma.Schema):
    structures = ma.fields.List(
        ma.fields.String(), required=True, validate=ma.validate.Length(min=1))
    include_embedding = ma.fields.Boolean(missing=False)


class ParticleSchema(ma.Schema):
    class Meta:
        unknown = ma.EXCLUDE

    structure = ma.fields.String(required=True)
    original_structure = ma.fields.String(required=False)
    # decoded = ma.fields.String(required=False)
    embedding = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(ma.fields.Float()), required=False)
    projection = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(
        ma.fields.Float()), required=False)
    properties = ma.fields.Raw(required=False)
    nearest_neighbors = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.Raw(), required=True)
    clusters = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.Raw(), required=True)


class ServerCollectionSchema(ma.Schema):
    data = ma.fields.Nested(ParticleSchema, many=True, required=True)
    projections = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.Raw(), required=True)


from .projection import ProjectionSchema, ProjectionArgsSchema, ProjectionWithModelsArgsSchema
from .pso import PSOArgsSchema
from .interpolation import InterpolationArgsSchema, InterpolatedParticleSchema, NeighborhoodSamplingArgsSchema
from .molecule import *
from .mmp import *
from .stoned_selfies import *

