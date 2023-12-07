import logging

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions

def update_props(mol):
    mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(mol)

def match_queries(mol, queries: list, react_atom_id: int):
    size = len(queries)
    result = np.zeros(size, dtype=int)
    for i, query in enumerate(queries):
        result[i] = match_reacting_atom(mol, query, react_atom_id)
    return result

def match_reacting_atom(mol, query, react_atom_id: int):
    try:
        for match_atoms in mol.GetSubstructMatches(query):
            # print('match_atoms:', match_atoms)
            if react_atom_id == match_atoms[0]:
                return True
        return False
    except Exception as e:
        logging.warning(f"{type(e).__name__}: {e} \n for query: {Chem.MolToSmarts(query)}\n for smiles: {Chem.MolToSmiles(mol)}")


def find_acid_carbon(acid, react_atoms: list):
    for i in react_atoms:
        if acid.GetAtomWithIdx(i).GetSymbol() == 'C':
            return i

def smarts_descriptor(rxn_smi, amine_queries: list, acid_queries: list):
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
    rxn.Initialize()
    amine, acid = rxn.GetReactants()
    update_props(amine)
    update_props(acid)
    react_atoms_ids = rxn.GetReactingAtoms()

    react_nitrogen = react_atoms_ids[0][0]
    react_carbon = find_acid_carbon(acid, react_atoms_ids[1])
    amine_descr = match_queries(amine, amine_queries, react_nitrogen)
    acid_descr = match_queries(acid, acid_queries, react_carbon)
    return amine_descr, acid_descr
