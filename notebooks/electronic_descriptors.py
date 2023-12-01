import sys
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import RDLogger

# Disable RDKit logging
RDLogger.DisableLog('rdApp.*')

import pandas as pd
import numpy as np

import dask
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

from linux_qm.src.util import check_amide_mapping, ind_rxn_map
from linux_qm.desc.electronic import orca_calculation, get_amine_electronic, get_acid_electronic


dask.config.set({"dataframe.convert-string": False})
ProgressBar().register()


def main():
    get_acid_electronic("[CH3:1][N+:2]([CH2:5][CH2:6][NH2:7])([CH3:4])[CH3:3].O[C:9]([CH2:11][CH2:12][O:13][CH2:14][CH2:15][O:16][CH2:17][CH2:18][O:19][CH2:20][CH2:21][O:22][CH2:23][CH2:24][N:25]=[N+:26]=[N-:27])=[O:10]>>[CH3:1][N+:2]([CH3:4])([CH3:3])[CH2:5][CH2:6][NH:7][C:9](=[O:10])[CH2:11][CH2:12][O:13][CH2:14][CH2:15][O:16][CH2:17][CH2:18][O:19][CH2:20][CH2:21][O:22][CH2:23][CH2:24][N:25]=[N+:26]=[N-:27]")

    df = pd.read_csv('../data/slv_amides/amide_training_set.csv')

    df['rxn_smi'] = df.amine_smi + '.' + df.acid_smi + '>>' + df.product_smi

    def ha_count(smi: str):
        mol = Chem.MolFromSmiles(smi)
        return mol.GetNumHeavyAtoms() if mol else None

    print('Filter by heavy atoms')
    ha = df.rxn_smi.str.split('>>').apply(lambda x: x[1]).apply(ha_count)
    df = df[ha < 40].copy()

    print('Mapping reactions')
    df.rxn_smi = df.rxn_smi.str.replace(r':\d+', '', regex=True)
    dds = dd.from_pandas(df.rxn_smi, npartitions=128)
    df.rxn_smi = dds.apply(ind_rxn_map, meta=dds).compute(scheduler='threads')
    df.dropna(subset=['rxn_smi'], inplace=True)

    # rxn_ok = df.rxn_smi.apply(check_amide_mapping)

    print('Cheking mappings')
    dds = dd.from_pandas(df.rxn_smi, npartitions=128)
    rxn_ok = dds.apply(check_amide_mapping, meta=dds).compute(scheduler='threads')
    df = df[rxn_ok].copy()

    print('Invalid mappings:', (~rxn_ok).sum())
    print('Reactions:', len(df))

    # print('Amine descriptor calculation:')
    # dds = dd.from_pandas(df.rxn_smi, npartitions=128)
    # amine_descr = dds.apply(get_amine_electronic,
    #                         # args=(query,),
    #                         meta=dds).compute(scheduler='processes')
    # df['amine_electronic_descr'] = amine_descr
    # print('Errors:', amine_descr.isna().sum())
    #
    # df.to_pickle('../data/slv_amides/amine_descriptor.pkl')

    print('Acid descriptor calculation:')
    dds = dd.from_pandas(df.rxn_smi, npartitions=128)
    acid_descr = dds.apply(get_acid_electronic,
                           # args=(query,),
                           meta=dds).compute(scheduler='processes')

    df['acid_electronic_descr'] = acid_descr
    print('Errors:', acid_descr.isna().sum())

    df.to_pickle('../data/slv_amides/acid_descriptor.pkl')

if __name__ == '__main__':
    main()