import marshmallow as ma


class MoleculeImageArgsSchema(ma.Schema):
    structure = ma.fields.String(required=True)
    substructure = ma.fields.String(required=False)


class MoleculesImageArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True, validate=ma.validate.Length(min=1))
    method = ma.fields.String(required=False, missing='mcs', validate=ma.validate.OneOf(['auto', 'single', 'murcko', 'mcs', 'similarity']))
    substructure = ma.fields.String(required=False)


class MoleculesSubstructureArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True, validate=ma.validate.Length(min=1))
    smarts = ma.fields.String(required=True, validate=ma.validate.Length(min=1))


class MoleculesSubstructureSchema(ma.Schema):
    validity = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.Boolean(), required=True)
    counts = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.Integer(), required=True)
    # substructures = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.List(ma.fields.Float()), required=True)


class MoleculesTanimotoArgsSchema(ma.Schema):
    structures = ma.fields.List(ma.fields.String(), required=True, validate=ma.validate.Length(min=1))
    reference = ma.fields.String(required=True, validate=ma.validate.Length(min=1))
    fingerprint = ma.fields.String(required=False, missing='daylight', validate=ma.validate.OneOf(['morgan', 'daylight', 'cddd']))


class MoleculesTanimotoSchema(ma.Schema):
    tanimoto = ma.fields.Dict(key=ma.fields.Str(), value=ma.fields.Boolean(), required=True)
