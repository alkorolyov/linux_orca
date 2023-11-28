import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions


if __name__ == '__main__':
    rxn_smi = sys.argv[1] if len(sys.argv) > 1 else 'N.COO>>NCO'

