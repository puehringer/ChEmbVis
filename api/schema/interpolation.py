import marshmallow as ma


class InterpolationArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True, validates=ma.validate.Length(min=2))
    maxSamples = ma.fields.Integer(strict=True, default=10, validate=ma.validate.Range(min=2))
