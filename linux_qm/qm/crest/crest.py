import re
import os
import sys
import shutil
import pickle
import subprocess
import logging
import argparse


from uuid import uuid4

from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%d-%m-%Y %I:%M:%S',
    encoding='utf-8',
    level=logging.INFO)

TMP_DIR = '.crest_tmp'

# def generate_constrain_input(constr_val_angles: list):
#     atoms_constr = []
#     for a in constr_val_angles:
#         atoms_constr.extend(mol.get_atoms_by_symbol(a))
#     cinp = "$constrain\n"
#     cinp += gen_angle_constraints(mol, atoms_constr)
#     cinp += "$end\n"
#
# def gen_angle_constraints(atoms: list):
#     """ Generate constraints for all angles where atom is the middle atom """
#     constr = ""
#     for atom_symbol in atoms:
#
#         neigbors = mol.get_connected_atoms(a)
#         for a1, a2 in combinations(neigbors, 2):
#             i1 = mol.get_atom_idx(a1) + 1
#             i2 = mol.get_atom_idx(a) + 1
#             i3 = mol.get_atom_idx(a2) + 1
#             angle = mol.get_angle(a1, a, a2) * 180 / pi
#             constr += f"  angle: {i1}, {i2}, {i3}, {angle:0.4f}\n"
#     return constr


def conformer_pipeline(smi: str, n_jobs: int):
    saved_work_dir = os.getcwd()

    tmp_path = _create_tmp_dir()
    os.chdir(tmp_path)

    mol = _load_smiles3D(smi)
    xyz_fname = _create_input_file(mol)

    _check_error(mol, xyz_fname, 'SMILES parsing')

    # pre-optimize
    xtb_optimize(
        xyz_fname,
        method="gff",
        n_jobs=n_jobs
    )

    _check_error(mol, 'xtbopt.xyz', 'geometry optimization')

    # generate conformers
    xyz_fname = 'xtbopt.xyz'
    conformer_gen(
        xyz_fname,
        method="gfnff",
        ewin=8.0,
        mdlen=5.0,
        n_jobs=n_jobs,
    )

    _check_error(mol, 'crest_conformers.xyz', 'conformer generation')

    # screen conformers
    xyz_fname = 'crest_conformers.xyz'
    conformer_screen(
        xyz_fname,
        method="gfn2",
        ewin=12.0,
        n_jobs=n_jobs,
        verbose=False)

    _check_error(mol, 'crest_ensemble.xyz', 'conformer screen')

    add_conformers(mol, 'crest_ensemble.xyz')

    # DEV

    original_mol = pickle.dumps(_load_smiles3D(smi))
    conf_mol = pickle.dumps(mol, protocol=4)
    conf_mol_5 = pickle.dumps(mol, protocol=5)
    print('original_pkl', len(original_mol))
    print('conf_pkl 4', len(conf_mol))
    print('conf_pkl 5', len(conf_mol_5))
    print('num conf:', mol.GetNumConformers())
    print('bytes per atom (conformers):', (len(conf_mol)- len(original_mol))/ (mol.GetNumConformers() + 1) / _load_smiles3D(smi).GetNumAtoms())
    # for conf in mol.GetConformers():
    #     print('3D:', conf.Is3D())
    #     print('energy:', conf.GetProp('energy'))
    #     print('id:', conf.GetId())


    os.chdir(saved_work_dir)
    shutil.rmtree(tmp_path)
    return mol

def _create_tmp_dir():
    uid = str(uuid4())
    tmp_path = f"{TMP_DIR}/{uid}"
    os.makedirs(tmp_path, exist_ok=True)
    return tmp_path

def _load_smiles3D(smi: str):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(mol)  # 3D coordinates
    mol.SetProp("smiles", Chem.CanonSmiles(smi))
    return mol

def _create_input_file(mol: Chem.rdchem.Mol):
    fname = 'input.xyz'
    xyz = Chem.MolToXYZBlock(mol)
    with open('input.xyz', 'w') as f:
        f.write(xyz)
    return fname

def _check_error(mol, filepath: str, stage: str):
    if os.path.exists(filepath):
        logging.debug(f'Success {stage}')
    else:
        raise ValueError(f"Failed {stage} for '{mol.GetProp('smiles')}'")



def xtb_optimize(
        xyz_name,
        method: str = "gff",
        crit: str = "normal",
        xtbinp: str = "",
        maxiter: int = 50,
        n_jobs: int = 8,
        verbose = False,
):
    """
    Attempt a geometry optimization with parameters specified
    """

    # command that will be used to execute xtb package
    _cmd = f"""xtb {xyz_name} --{method} --opt {crit} --cycles {maxiter} {"--input param.inp" if xtbinp else ""} -P {n_jobs}"""
    # print(_cmd)
    subprocess.run(
        _cmd.split(' '),
        stdout = None if verbose else subprocess.DEVNULL,
    )

def conformer_gen(
        xyz_name,
        scratch = False,
        method: str = "gfnff",
        ewin: float = 6.0,
        mdlen: float = 15.0,
        mddump: float = 250.0,
        vbdump: float = 1.0,
        chk_topo = False,
        constr_val_angles: list = [],
        force_const: float = 0.05,
        n_jobs: int = 8,
        verbose = False,
):
    _cmd = f"""crest {xyz_name} {f'--scratch' if scratch else ''} -{method} -quick -ewin {ewin:0.4f} -mdlen {mdlen:0.4f} -mddump {mddump:0.4f} -vbdump {vbdump:0.4f} -T {n_jobs}"""

    if not chk_topo:
        _cmd += " --noreftopo"

    # TODO implement valent angles constrains
    # if constr_val_angles:
    #     _cmd += f" -cinp {name}_constraint.inp -fc {force_const:0.4f}"

    # print(_cmd)
    # print(_cmd)
    subprocess.run(
        _cmd.split(' '),
        stdout = None if verbose else subprocess.DEVNULL,
    )

def conformer_screen(
        xyz_name: str,
        method: str = "gfn2",
        ewin: float = 6.0,
        n_jobs: int = 8,
        verbose = False,
):
    """
    Any conformer ensemble present in a molecule is reoptimized with the selected method,
    then pruned and sorted by energies. Useful in eliminating redundancies from deficiencies of GFN-FF, for instance.
    """

    _cmd = f"""crest -screen {xyz_name} -{method} -ewin {ewin} -T {n_jobs} """

    subprocess.run(
        _cmd.split(' '),
        stdout = None if verbose else subprocess.DEVNULL,
    )



def add_conformers(mol, ensemble_path: str):
    # remove existing conf
    mol.RemoveConformer(0)

    conformers = _read_conformers(ensemble_path)

    for conf in conformers:
        energy, coords = _parse_xyz(conf)
        _add_conformer_to_molecule(mol, energy, coords)

def _read_conformers(ensemble_path):
    conformers = []
    with open(ensemble_path, 'r') as file:
        lines = file.readlines()

    conformer_str = ''
    for line in lines:
        if re.match(r'^\s*\d+', line):  # beginning of new conformer
            if conformer_str:
                conformers.append(conformer_str)
            conformer_str = line    # new conformer string
        else:
            conformer_str += line

    if conformer_str:   # add last conformer
        conformers.append(conformer_str)

    return conformers

def _add_conformer_to_molecule(molecule, energy, coordinates):
    conformer = Chem.Conformer()
    conformer.SetProp('energy', str(energy))
    for idx, (x, y, z) in enumerate(coordinates):
        conformer.SetAtomPosition(idx, (x, y, z))
    molecule.AddConformer(conformer, assignId=True)

def _parse_xyz(xyz: str):
    lines = xyz.split('\n')

    # Extract atom count and coordinates from XYZ block
    atom_count = int(lines[0].strip())
    energy = float(lines[1].strip())
    coordinates = [line.split() for line in lines[2:2 + atom_count]]

    # Convert coordinates to RDKit format
    atom_symbols = [coord[0] for coord in coordinates]
    xyz_coordinates = [(float(coord[1]), float(coord[2]), float(coord[3])) for coord in coordinates]

    return energy, xyz_coordinates


def main():
    parser = argparse.ArgumentParser(description='Generate conformer ensemble using CREST')
    #
    # Define command-line arguments
    parser.add_argument('input_smiles', help='Input structure SMILES')
    parser.add_argument('--n_jobs', '-n', default=8, type=int, help='Number of parallel threads to use')
    parser.add_argument('--verbose', '-v', default=False, action='store_true', help='Enable verbose mode')

    # Parse the command-line arguments
    args = parser.parse_args()

    smiles = args.input_smiles
    n_jobs = args.n_jobs
    verbose = args.verbose

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    mol = conformer_pipeline(smiles, n_jobs)

    return mol


if __name__ == '__main__':
    # smi = 'COC1=CC=CC([C@](O2)(CN3C=CN=C3)OC[C@H]2COC4=CC=CC=C4)=C1'
    smi = 'COC1CN(C(OC)COC)C1'
    # smi = 'COC1CC(NC)C1'
    # smi = 'COC1CNC1'
    # smi = 'O'

    # DEV
    sys.argv[1:] = [smi,'-n 8', '-v']

    main()

