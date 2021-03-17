import marshmallow as ma

from .projection import ProjectionSchema, ProjectionArgsSchema
from .pso import PSOArgsSchema
from .interpolation import InterpolationArgsSchema
from .molecule import *
from .mmp import *


class EmbeddingArgsSchema(ma.Schema):
    structures = ma.fields.List(
        ma.fields.String(), required=True, validate=ma.validate.Length(min=1))


class ParticleSchema(ma.Schema):
    structure = ma.fields.String(required=True)
    # decoded = ma.fields.String(required=False)
    embedding = ma.fields.List(ma.fields.Float(), required=True)
    projection = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(
        ma.fields.Float()), required=False)
    properties = ma.fields.Raw(required=False)
