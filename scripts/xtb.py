import os
import shutil
import logging
import subprocess
from time import time
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from linux_qm.src.util import load_smiles3D, _create_tmp_dir, draw3Dconfs
from linux_qm.qm.driver import set_positions
from linux_qm.qm.crest.crest import conformer_pipeline

smi = 'COC1CN(C)C1'
# smi = 'CO'


mol = load_smiles3D(smi)
EmbedMultipleConfs(mol, numConfs=3)

# mol = conformer_pipeline(smi)

# Get the conformers

conformers = mol.GetConformers()

# Print the number of conformers generated
print(f"Number of conformers generated: {len(conformers)}")
# print("Energies:", [c.GetProp('energy') for c in conformers])

from linux_qm.qm.orca.orca import OrcaDriver
orca = OrcaDriver()

logging.getLogger().setLevel(logging.DEBUG)
# orca.options['solvent'] = 'THF'

orca.options['method'] = 'HF-3c'
orca.options['n_jobs'] = 1
conf = deepcopy(mol).GetConformer()
orca.single_point(conf)


# conf = deepcopy(mol).GetConformer()
# orca.options['method'] = 'M062X cc-pVDZ NMR'
# orca.geometry_optimization(conf)

# conf = deepcopy(mol).GetConformer()
# orca.options['method'] = 'M062X aug-cc-pVTZ NMR'
# orca.geometry_optimization(conf)

# conf = deepcopy(mol).GetConformer()
# orca.options['method'] = 'M062X aug-cc-pVTZ NMR'
# orca.options['n_jobs'] = 48
# orca.geometry_optimization(conf)

# orca.options['n_jobs'] = 16