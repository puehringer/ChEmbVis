import marshmallow as ma


class PSOArgsSchema(ma.Schema):
    structure = ma.fields.String(
        required=True, validate=ma.validate.Length(min=2))
    num_swarms = ma.fields.Integer(
        required=False, missing=3, validate=ma.validate.Range(min=1))
    num_part = ma.fields.Integer(
        required=False, missing=30, validate=ma.validate.Range(min=1))
    iterations = ma.fields.Integer(
        required=False, missing=50, validate=ma.validate.Range(min=1))
    v_min = ma.fields.Float(required=False, missing=None)
    v_max = ma.fields.Float(required=False, missing=None)
    inertia_weight = ma.fields.Float(required=False, missing=None)
    phi1 = ma.fields.Float(required=False, missing=None)
    phi2 = ma.fields.Float(required=False, missing=None)
    phi3 = ma.fields.Float(required=False, missing=None)
    # TODO: Create model
    objectives = ma.fields.List(ma.fields.Raw(), required=True, validate=ma.validate.Length(min=1))
