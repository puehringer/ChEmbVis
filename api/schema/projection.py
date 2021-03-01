import marshmallow as ma


class ProjectionArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True)
    embedding = ma.fields.List(ma.fields.List(ma.fields.Float(), validate=ma.validate.Length(min=2)), required=True, validate=ma.validate.Length(min=2))
    additional_structures = ma.fields.List(ma.fields.String(), required=False)
    additional_embedding = ma.fields.List(ma.fields.List(ma.fields.Float(), validate=ma.validate.Length(min=2)), required=False, validate=ma.validate.Length(min=2))


class ProjectionSchema(ma.Schema):
    projections = ma.fields.List(ma.fields.String(), required=True)
    projection = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(ma.fields.Float()), required=True)
    additional = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(ma.fields.Float()), required=False)
