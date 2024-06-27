import sys
import logging
import re

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

import pandas as pd
import numpy as np

from linux_qm.src.render import indigo, draw_reacting_mapnums
from linux_qm.src.util import load_smiles3D, ind_rxn_map
from linux_qm.qm.orca.orca import OrcaDriver

def orca_calculation(conf):
    orca = OrcaDriver()

    orca.options['n_jobs'] = 1

    orca.options['method'] = 'XTB2'
    orca.geometry_optimization(conf)

    orca.options['method'] = 'HF-3c'
    # orca.options['method'] = 'BP86 def2-SVP def2/J D3BJ RIJCOSX'
    # orca.options['method'] = 'BP86 def2-TZVP def2/J D3BJ RIJCOSX'
    data = orca.single_point(conf)
    # data = orca.geometry_optimization(conf, calc_npa=False)

    # data = orca.geometry_optimization(conf, calc_npa=True)
    return data

def parse_qm_results(data, atom_ids):
    mull_charges = data.atomcharges['mulliken'][atom_ids]
    low_charges = data.atomcharges['lowdin'][atom_ids]
    # npa_charges = data.atomcharges['npa'][amine_raids]

    homo = data.moenergies[0][data.homos[0]]
    lumo = data.moenergies[0][data.homos[0] + 1]

    logging.debug(f"Mulliken charges: {mull_charges}")
    logging.debug(f"Lowdig charges: {low_charges}")
    logging.debug(f"HOMO LUMO: {homo} {lumo}")

    descriptor = np.hstack([mull_charges, low_charges, homo, lumo]).round(6)
    return descriptor

def get_amine_atoms(mol, reacting_nitrogen_idx: int):
    amine_aids = [reacting_nitrogen_idx]
    atom = mol.GetAtomWithIdx(reacting_nitrogen_idx)
    for a in atom.GetNeighbors():
        if a.GetSymbol() == 'C':
            amine_aids.append(a.GetIdx())
    return amine_aids

def get_amine_electronic(rxn_smi, automap=False):
    try:
        if automap:
            rxn_smi = ind_rxn_map(rxn_smi)

        # load rxn
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
        rxn.Initialize()

        # get reactants
        amine, acid = rxn.GetReactants()
        amine_atoms = get_amine_atoms(amine, rxn.GetReactingAtoms()[0][0])

        mol = load_smiles3D(Chem.MolToSmiles(amine), opt=True)

        # qm calculation
        data = orca_calculation(mol.GetConformer())

        for i in amine_atoms:
            logging.debug(amine.GetAtomWithIdx(i).GetSymbol())

        return parse_qm_results(data, amine_atoms)

    except Exception as e:
        logging.warning(f"{type(e).__name__}: {e} for rxn_smi: {rxn_smi}")

def get_acid_atoms(mol, react_acid_ids: list):
    # Order should be
    # [C, =O, -OH, C bounded to C=O]
    result = []
    for i in react_acid_ids:
        if mol.GetAtomWithIdx(i).GetSymbol() == 'C':
            result.append(i)
    for i in react_acid_ids:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
            result.append(i)
    for i in react_acid_ids:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 0:
            result.append(i)

    atom = mol.GetAtomWithIdx(result[0])
    for nei in atom.GetNeighbors():
        if nei.GetSymbol() == 'C':
            result.append(nei.GetIdx())

    return result

def get_acid_electronic(rxn_smi, automap=False):
    try:
        if automap:
            rxn_smi = ind_rxn_map(rxn_smi)

        # load rxn
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
        rxn.Initialize()

        # get reactants
        amine, acid = rxn.GetReactants()
        Chem.SanitizeMol(acid)
        Chem.SanitizeMol(amine)

        acid_atoms = get_acid_atoms(acid, rxn.GetReactingAtoms()[1])
        logging.debug(f'Acid atom ids: {acid_atoms}')

        logging.info(f"Heavy Atom Count: {acid.GetNumHeavyAtoms()}")

        mol = load_smiles3D(Chem.MolToSmiles(acid), opt=True)

        # qm calculation
        data = orca_calculation(mol.GetConformer())

        return parse_qm_results(data, acid_atoms)

    except Exception as e:
        logging.warning(f"{type(e).__name__}: {e} for rxn_smi: {rxn_smi}")




