import sys
import shutil

from rdkit import Chem
from rdkit.Chem import AllChem

from linux_qm.src.util import _load_smiles3D

def main():
    pass

if __name__ == '__main__':
    smi = 'COC1CNC1'

    mol = _load_smiles3D(smi)


    sys.argv[1:] = ['']
    main()

