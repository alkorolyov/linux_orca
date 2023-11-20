from typing import List
from rdkit import Chem

class Conformer:
    energy: float
    xyz: str
    qm_method: str
    descriptors: list


class Molecule:
    name: str
    smiles: str
    rdmol: Chem.rdchem.Mol
    conformers: List[Conformer]


def parse_crest_conformers(original_mol, conformers_fpath: str):
    with open(conformers_fpath, 'r') as f:
        xyz_buffer = f.readlines()








