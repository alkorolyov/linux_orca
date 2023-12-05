import sys
import os
import shutil
import logging
import subprocess
from copy import deepcopy
from time import time

from rdkit import Chem
from rdkit.Chem import AllChem

from linux_qm.src.util import load_smiles3D
from linux_qm.qm.crest.crest import conformer_pipeline
from linux_qm.qm.orca.orca import OrcaDriver

logging.getLogger().setLevel(logging.DEBUG)

if __name__=='__main__':

    if len(sys.argv) == 1:
        sys.argv[1:] = [1]

    # smi = 'COC1=CC=CC([C@](O2)(CN3C=CN=C3)OC[C@H]2COC4=CC=CC=C4)=C1'
    smi = 'COC1CN(C(OC)COC)C1'
    # smi = 'CCCOC1CN(C)C1'
    # smi = 'COC1CN(C)C1'
    # smi = 'CO'

    n_jobs = int(sys.argv[1])

    # mol = load_smiles3D(smi)

    mol = conformer_pipeline(smi, n_jobs=n_jobs)

    #
    # conf = mol.GetConformer()
    #
    # options = {
    #     'method': 'B3LYP 6-31G',
    #     'solvent': None,
    #     'n_jobs': n_jobs,
    # }
    # orca = OrcaDriver(options=options)
    # orca.geometry_optimization(conf)