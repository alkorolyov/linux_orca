import sys
import logging
import re

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions


import pandas as pd
import numpy as np

from linux_qm.src.render import indigo, draw_reacting_mapnums
from linux_qm.src.util import load_smiles3D
from linux_qm.qm.orca.orca import OrcaDriver

from rxnmapper import BatchedMapper

def rxn_map_arr(rxn_smi_array: pd.Series, batch_size=8):
    rxn_mapper = BatchedMapper(batch_size=8)
    rxn_smi_clear = rxn_smi_array.str.replace(r':\d+','', regex=True)
    mapped_rxn = list(rxn_mapper.map_reactions(rxn_smi_clear))
    return mapped_rxn


def rxn_map(rxn_smi):
    rxn_mapper = BatchedMapper(batch_size=1)
    clear_smi = re.sub(r':\d+','',  rxn_smi)
    mapped_rxn = list(rxn_mapper.map_reactions([clear_smi]))[0]
    return mapped_rxn


def get_amine_atoms(mol, reacting_nitrogen_idx: int):
    amine_aids = [reacting_nitrogen_idx]
    atom = mol.GetAtomWithIdx(reacting_nitrogen_idx)
    for a in atom.GetNeighbors():
        if a.GetSymbol() == 'C':
            amine_aids.append(a.GetIdx())
    return amine_aids

def orca_calculation(conf):
    orca = OrcaDriver()

    orca.options['n_jobs'] = 1

    orca.options['method'] = 'XTB2'
    orca.geometry_optimization(conf)

    orca.options['method'] = 'HF-3c'
    # orca.options['method'] = 'BP86 def2-SVP def2/J D3BJ RIJCOSX'
    data = orca.single_point(conf, calc_npa=False)

    # data = orca.geometry_optimization(conf, calc_npa=True)
    return data

def gen_amine_electronic(rxn_smi):
    # try rxn_mapper:
    rxn_smi = rxn_map(rxn_smi)

    # load rxn
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
    rxn.Initialize()

    # get reactants
    amine, acid = rxn.GetReactants()
    amine_raids = get_amine_atoms(amine, rxn.GetReactingAtoms()[0][0])
    logging.debug(f'Amine atom ids: {amine_raids}')

    logging.info(f"Heavy Atom Count: {amine.GetNumHeavyAtoms()}")

    mol = load_smiles3D(Chem.MolToSmiles(amine), opt=True)

    # qm calculation
    data = orca_calculation(mol.GetConformer())

    charges = np.hstack([
        data.atomcharges['mulliken'][amine_raids],
        data.atomcharges['lowdin'][amine_raids],
        # data.atomcharges['npa'][amine_raids],
    ])
    homo = data.moenergies[0][data.homos[0]]
    lumo = data.moenergies[0][data.homos[0] + 1]

    return np.hstack([charges, homo, lumo]).round(6)