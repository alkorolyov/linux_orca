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
    try:
        # load rxn
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
        rxn.Initialize()

        # get reactants
        amine, acid = rxn.GetReactants()
        amine_raids = get_amine_atoms(amine, rxn.GetReactingAtoms()[0][0])

        mol = load_smiles3D(Chem.MolToSmiles(amine), opt=True)

        # qm calculation
        data = orca_calculation(mol.GetConformer())

        for i in amine_raids:
            logging.debug(amine.GetAtomWithIdx(i).GetSymbol())

        mull_charges = data.atomcharges['mulliken'][amine_raids]
        low_charges = data.atomcharges['lowdin'][amine_raids]
        # npa_charges = data.atomcharges['npa'][amine_raids]

        homo = data.moenergies[0][data.homos[0]]
        lumo = data.moenergies[0][data.homos[0] + 1]

        logging.debug(f"Mulliken charges: {mull_charges}")
        logging.debug(f"Lowdig charges: {low_charges}")
        logging.debug(f"HOMO LUMO: {homo} {lumo}")

        descriptor = np.hstack([mull_charges, low_charges, homo, lumo]).round(6)
        return descriptor

    except Exception as e:
        logging.warning(f"{type(e).__name__}: {e} for rxn_smi: {rxn_smi}")

