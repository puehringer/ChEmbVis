import marshmallow as ma
from . import ParticleSchema

class InterpolationArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True, validate=ma.validate.Length(min=2))
    maxSamples = ma.fields.Integer(strict=True, missing=10, validate=ma.validate.Range(min=2))


class InterpolatedParticleSchema(ParticleSchema):
    image = ma.fields.String(required=False)
    image1 = ma.fields.String(required=False)