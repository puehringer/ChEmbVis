import marshmallow as ma


class InterpolationArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True, validate=ma.validate.Length(min=2))
    maxSamples = ma.fields.Integer(strict=True, missing=10, validate=ma.validate.Range(min=2))
