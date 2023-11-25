import sys
import os
import shutil
import logging
import subprocess
from copy import deepcopy
from time import time

from linux_qm.qm.crest.crest import conformer_pipeline
from linux_qm.qm.orca.orca import OrcaDriver

# logging.getLogger().setLevel(logging.DEBUG)

if __name__=='__main__':
    # smi = 'COC1CN(C)C1'
    smi = 'CO'
    n_jobs = sys.argv[1]
    mol = conformer_pipeline(smi, n_jobs=n_jobs)
    conf = mol.GetConformer()

    options = {
        'method': 'M062X 6-31G NMR',
        'solvent': 'THF',
        'n_jobs': n_jobs,
    }
    orca = OrcaDriver(options=options)
    orca.geometry_optimization(conf)