import os
from uuid import uuid4

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem


def load_smiles3D(smi: str, opt=False, num_conf=1):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)  # 3D coordinates
    if num_conf > 1:
        AllChem.EmbedMultipleConfs(mol, num_conf)
    if opt:
        AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp("smiles", Chem.CanonSmiles(smi))
    return mol


def SetPositions(conf, atom_positions):
    """
    Example usage:
    atom_positions = [
        [2.225, -0.136, -0.399],
        [1.158, -0.319, 0.424],
        [-0.050, 0.113, 0.042],
        [ 0.178,-0.956, 0.329],
    ]
    conf = mol.GetConformer()
    SetPositions(conf, atom_positions)

    conf: RDKit conformer
    atom_positions: 2D array of x, y, z coordinates of each atom
    """

    for i in range(conf.GetNumAtoms()):
        xyz = atom_positions[i]
        conf.SetAtomPosition(i, xyz)

def draw3Dconfs(mol, autoalign=True, confIds=None, size=(600, 400)):
    if autoalign:
        AllChem.AlignMolConformers(mol)

    print('num conformers', mol.GetNumConformers())

    # Visualize using Py3Dmol
    viewer = py3Dmol.view(width=size[0], height=size[1])

    if not confIds:
        # all conformers
        confIds = [conf.GetId() for conf in mol.GetConformers()]

    for conf in mol.GetConformers():
        if conf.GetId() in confIds:
            mol_block = Chem.MolToMolBlock(mol, confId=conf.GetId())
            viewer.addModel(mol_block, "mol")

    viewer.setStyle({"stick": {}})
    # viewer.setBackgroundColor("black")
    viewer.zoomTo()
    viewer.show()



def _create_tmp_dir(tmp_root):
    uid = str(uuid4())
    tmp_path = f"{tmp_root}/{uid}"
    os.makedirs(tmp_path, exist_ok=True)
    return tmp_path
