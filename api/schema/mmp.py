import marshmallow as ma


class MMPArgsSchema(ma.Schema):
    structure = ma.fields.String(required=True)
    min_variable_size = ma.fields.Int(missing=0)
    max_variable_size = ma.fields.Int(missing=9999)
    min_constant_size = ma.fields.Int(missing=0)
    min_radius = ma.fields.Int(missing=0)
    min_pairs = ma.fields.Int(missing=0)
    substructure = ma.fields.String(missing='')


class MMPSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True)

