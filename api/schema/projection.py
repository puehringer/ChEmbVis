import marshmallow as ma
from . import ParticleSchema


class ProjectionArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True)
    embedding = ma.fields.List(ma.fields.List(ma.fields.Float(), validate=ma.validate.Length(min=2)), required=True, validate=ma.validate.Length(min=2))
    additional_structures = ma.fields.List(ma.fields.String(), required=False)
    additional_embedding = ma.fields.List(ma.fields.List(ma.fields.Float(), validate=ma.validate.Length(min=2)), required=False, validate=ma.validate.Length(min=2))


class ProjectionWithModelsParticlesSchema(ma.Schema):
    class Meta:
        unknown = ma.EXCLUDE
    
    structure = ma.fields.String(required=True)
    embedding = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(ma.fields.Float()), required=True)


class ProjectionWithModelsArgsSchema(ma.Schema):
    class Meta:
        unknown = ma.EXCLUDE

    particles = ma.fields.Nested(ProjectionWithModelsParticlesSchema, many=True, required=True)
    models = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.Str(), required=False)


class ProjectionSchema(ma.Schema):
    projections = ma.fields.Dict(key=ma.fields.String(), value=ma.fields.Raw(), required=True)
    projection = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(ma.fields.Float()), required=True)
    additional = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(ma.fields.Float()), required=False)
