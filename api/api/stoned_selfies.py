import time
import selfies
import rdkit
import random
import numpy as np
import random
from rdkit import Chem
from selfies import encoder, decoder
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem import Mol
from rdkit.Chem.AtomPairs.Sheridan import GetBPFingerprint, GetBTFingerprint
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem import Draw

from rdkit.Chem import MolToSmiles as mol2smi
from rdkit import RDLogger
from flask import jsonify
from ..utils import cached, mol, logger, catch_time
from ..constants import blp
from ..schema import StonedSelfiesArgsSchema, StonedSelfiesSchema

RDLogger.DisableLog('rdApp.*')


def randomize_smiles(mol):
    '''Returns a random (dearomatized) SMILES given an rdkit mol object of a molecule.
    Parameters:
    mol (rdkit.Chem.rdchem.Mol) :  RdKit mol object (None if invalid smile string smi)

    Returns:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object  (None if invalid smile string smi)
    '''
    if not mol:
        return None

    Chem.Kekulize(mol)
    return rdkit.Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False,  kekuleSmiles=True)


def sanitize_smiles(smi):
    '''Return a canonical smile representation of smi

    Parameters:
    smi (string) : smile string to be canonicalized 

    Returns:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object                          (None if invalid smile string smi)
    smi_canon (string)          : Canonicalized smile representation of smi (None if invalid smile string smi)
    conversion_successful (bool): True/False to indicate if conversion was  successful 
    '''
    try:
        mol = smi2mol(smi, sanitize=True)
        smi_canon = mol2smi(mol, isomericSmiles=False, canonical=True)
        return (mol, smi_canon, True)
    except:
        return (None, None, False)


def get_selfie_chars(selfie):
    '''Obtain a list of all selfie characters in string selfie

    Parameters: 
    selfie (string) : A selfie string - representing a molecule 

    Example: 
    >>> get_selfie_chars('[C][=C][C][=C][C][=C][Ring1][Branch1_1]')
    ['[C]', '[=C]', '[C]', '[=C]', '[C]', '[=C]', '[Ring1]', '[Branch1_1]']

    Returns:
    chars_selfie: list of selfie characters present in molecule selfie
    '''
    chars_selfie = []  # A list of all SELFIE sybols from string selfie
    while selfie != '':
        chars_selfie.append(selfie[selfie.find('['): selfie.find(']')+1])
        selfie = selfie[selfie.find(']')+1:]
    return chars_selfie


class _FingerprintCalculator:
    ''' Calculate the fingerprint for a molecule, given the fingerprint type
    Parameters: 
        mol (rdkit.Chem.rdchem.Mol) : RdKit mol object (None if invalid smile string smi)
        fp_type (string)            :Fingerprint type  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6)  
    Returns:
        RDKit fingerprint object
    '''

    def get_fingerprint(self, mol: Mol, fp_type: str):
        method_name = 'get_' + fp_type
        method = getattr(self, method_name)
        if method is None:
            raise Exception(f'{fp_type} is not a supported fingerprint type.')
        return method(mol)

    def get_AP(self, mol: Mol):
        return AllChem.GetAtomPairFingerprint(mol, maxLength=10)

    def get_PHCO(self, mol: Mol):
        return Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)

    def get_BPF(self, mol: Mol):
        return GetBPFingerprint(mol)

    def get_BTF(self, mol: Mol):
        return GetBTFingerprint(mol)

    def get_PATH(self, mol: Mol):
        return AllChem.RDKFingerprint(mol)

    def get_ECFP4(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 2)

    def get_ECFP6(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 3)

    def get_FCFP4(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 2, useFeatures=True)

    def get_FCFP6(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 3, useFeatures=True)


def get_fingerprint(mol: Mol, fp_type: str):
    ''' Fingerprint getter method. Fingerprint is returned after using object of 
        class '_FingerprintCalculator'

    Parameters: 
        mol (rdkit.Chem.rdchem.Mol) : RdKit mol object (None if invalid smile string smi)
        fp_type (string)            :Fingerprint type  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6)  
    Returns:
        RDKit fingerprint object

    '''
    return _FingerprintCalculator().get_fingerprint(mol=mol, fp_type=fp_type)


def mutate_selfie(selfie: str, max_molecules_len: int, substructure: str = None, write_fail_cases: bool = False):
    '''Return a mutated selfie string (only one mutation on slefie is performed)

    Mutations are done until a valid molecule is obtained 
    Rules of mutation: With a 33.3% propbabily, either: 
        1. Add a random SELFIE character in the string
        2. Replace a random SELFIE character with another
        3. Delete a random character

    Parameters:
    selfie            (string)  : SELFIE string to be mutated 
    max_molecules_len (int)     : Mutations of SELFIE string are allowed up to this length
    write_fail_cases  (bool)    : If true, failed mutations are recorded in "selfie_failure_cases.txt"

    Returns:
    selfie_mutated    (string)  : Mutated SELFIE string
    smiles_canon      (string)  : canonical smile of mutated SELFIE string
    '''
    valid = False
    fail_counter = 0
    chars_selfie = get_selfie_chars(selfie)

    while not valid:
        fail_counter += 1

        # 34 SELFIE characters
        alphabet = list(selfies.get_semantic_robust_alphabet())

        choice_ls = [1, 2, 3]  # 1=Insert; 2=Replace; 3=Delete
        random_choice = np.random.choice(choice_ls, 1)[0]

        # Insert a character in a Random Location
        if random_choice == 1:
            random_index = np.random.randint(len(chars_selfie)+1)
            random_character = np.random.choice(alphabet, size=1)[0]

            selfie_mutated_chars = chars_selfie[:random_index] + \
                [random_character] + chars_selfie[random_index:]

        # Replace a random character
        elif random_choice == 2:
            random_index = np.random.randint(len(chars_selfie))
            random_character = np.random.choice(alphabet, size=1)[0]
            if random_index == 0:
                selfie_mutated_chars = [
                    random_character] + chars_selfie[random_index+1:]
            else:
                selfie_mutated_chars = chars_selfie[:random_index] + [
                    random_character] + chars_selfie[random_index+1:]

        # Delete a random character
        elif random_choice == 3:
            random_index = np.random.randint(len(chars_selfie))
            if random_index == 0:
                selfie_mutated_chars = chars_selfie[random_index+1:]
            else:
                selfie_mutated_chars = chars_selfie[:random_index] + \
                    chars_selfie[random_index+1:]

        else:
            raise Exception('Invalid Operation trying to be performed')

        selfie_mutated = "".join(x for x in selfie_mutated_chars)
        sf = "".join(x for x in chars_selfie)

        try:
            smiles = decoder(selfie_mutated)
            mol, smiles_canon, done = sanitize_smiles(smiles)
            if len(selfie_mutated_chars) > max_molecules_len or smiles_canon == "" or (substructure and not mol.HasSubstructMatch(rdkit.Chem.MolFromSmarts(substructure))):
                done = False
            if done:
                valid = True
            else:
                valid = False
        except:
            valid = False
            if fail_counter > 1 and write_fail_cases == True:
                f = open("selfie_failure_cases.txt", "a+")
                f.write('Tried to mutate SELFIE: '+str(sf) +
                        ' To Obtain: '+str(selfie_mutated) + '\n')
                f.close()

    return (selfie_mutated, smiles_canon)


def get_mutated_SELFIES(selfies_ls, num_mutations, substructure: str = None):
    ''' Mutate all the SELFIES in 'selfies_ls' 'num_mutations' number of times. 

    Parameters:
    selfies_ls   (list)  : A list of SELFIES 
    num_mutations (int)  : number of mutations to perform on each SELFIES within 'selfies_ls'
    substructure (str)  : Optional substructure to be matched.

    Returns:
    selfies_ls   (list)  : A list of mutated SELFIES

    '''
    for _ in range(num_mutations):
        selfie_ls_mut_ls = []
        for str_ in selfies_ls:

            str_chars = get_selfie_chars(str_)
            max_molecules_len = len(str_chars) + num_mutations

            selfie_mutated, _ = mutate_selfie(str_, max_molecules_len, substructure)
            selfie_ls_mut_ls.append(selfie_mutated)

        selfies_ls = selfie_ls_mut_ls.copy()
    return selfies_ls


def get_fp_scores(smiles_back, target_smi, fp_type):
    '''Calculate the Tanimoto fingerprint (using fp_type fingerint) similarity between a list 
       of SMILES and a known target structure (target_smi). 

    Parameters:
    smiles_back   (list) : A list of valid SMILES strings 
    target_smi (string)  : A valid SMILES string. Each smile in 'smiles_back' will be compared to this stucture
    fp_type (string)     : Type of fingerprint  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6) 

    Returns: 
    smiles_back_scores (list of floats) : List of fingerprint similarities
    '''
    smiles_back_scores = []
    target = Chem.MolFromSmiles(target_smi)

    fp_target = get_fingerprint(target, fp_type)

    for item in smiles_back:
        mol = Chem.MolFromSmiles(item)
        fp_mol = get_fingerprint(mol, fp_type)
        score = TanimotoSimilarity(fp_mol, fp_target)
        smiles_back_scores.append(score)
    return smiles_back_scores


def get_stoned_selfies(structure: str, substructure: str = None, random_samples: int = 1000, max_mutations: int = 5):
    # num_random_samples = 50000 # For a more exhaustive search!

    with catch_time('Randomizing molecules (in SELFIES)'):
        num_mutation_ls = list(range(max_mutations))

        mol = Chem.MolFromSmiles(structure)
        if mol == None:
            raise Exception('Invalid starting structure encountered')

        randomized_smile_orderings = [
            randomize_smiles(mol) for _ in range(random_samples)]

        # Convert all the molecules to SELFIES
        selfies_ls = [encoder(x) for x in randomized_smile_orderings]

    with catch_time('Mutation obtainment time (back to smiles)'):
        all_smiles_collect = []
        all_smiles_collect_broken = []

        for num_mutations in num_mutation_ls:
            # Mutate the SELFIES:
            selfies_mut = get_mutated_SELFIES(
                selfies_ls.copy(), num_mutations=num_mutations, substructure=substructure)

            # Convert back to SMILES:
            smiles_back = [decoder(x) for x in selfies_mut]
            all_smiles_collect = all_smiles_collect + smiles_back
            all_smiles_collect_broken.append(smiles_back)

    with catch_time('Unique mutated structure obtainment time'):
        # Work on:  all_smiles_collect
        canon_smi_ls = []
        for item in all_smiles_collect:
            mol, smi_canon, did_convert = sanitize_smiles(item)
            if mol == None or smi_canon == '' or did_convert == False:
                raise Exception('Invalid smile string found')
            canon_smi_ls.append(smi_canon)
        canon_smi_ls = list(set(canon_smi_ls))

    return canon_smi_ls


@blp.route('/stoned_selfies/', methods=['POST'])
@blp.arguments(StonedSelfiesArgsSchema, location='json')
@blp.response(200, StonedSelfiesSchema)
@cached
def stoned_selfies_route(args):
    structure = args.get('structure')
    substructure = args.get('substructure')
    random_samples = args.get('random_samples')
    max_mutations = args.get('max_mutations')

    canon_smi_ls = get_stoned_selfies(structure, random_samples=random_samples, max_mutations=max_mutations, substructure=substructure)

    with catch_time('Fingerprint calculation time'):
        canon_smi_ls_scores = get_fp_scores(canon_smi_ls, target_smi=structure, fp_type='ECFP4')

        filtered_smiles = [smi for smi, similarity in zip(canon_smi_ls, canon_smi_ls_scores) if similarity >= 0.0]

    # fingerprint = args.get('fingerprint')
    # scale = args.get('scale')
    return {"structures": filtered_smiles}
