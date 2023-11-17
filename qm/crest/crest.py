import os
import sys
import shutil
import subprocess

from time import time

# import molli as ml
# from indigo import Indigo
# from openbabel import openbabel as ob
from rdkit.Chem import AllChem as Chem


TMP_DIR = '.crest_tmp'

# obmol = ob.OBMol()
# obconv = ob.OBConversion()
# obconv.SetInAndOutFormats("smi", "mol2")
# obconv.AddOption("h", ob.OBConversion.GENOPTIONS)
# gen3d = ob.OBOp.FindType("gen3D")
#
# indigo = Indigo()

def smi2xyz(smi):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    Chem.EmbedMolecule(mol)
    return Chem.MolToXYZBlock(mol)

# def smi_to_mol2(smi):
#     obconv.ReadString(obmol, smi)
#     gen3d.Do(obmol, "--fastest")
#     return obconv.WriteString(obmol)
#
# def mol2xyz(mol):
#     xyz = f'{mol.countAtoms()}\n\n'
#     for atom in mol.iterateAtoms():
#         # print(atom.symbol(), atom.index())
#         x, y, z = atom.xyz()
#         xyz += f"{atom.symbol()}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
#     return xyz

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
        n_jobs: int = 12,
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
        n_jobs: int = 12,
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


def xtb_optimize(
        xyz_name,
        method: str = "gff",
        crit: str = "normal",
        xtbinp: str = "",
        maxiter: int = 50,
        n_jobs: int = 1,
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


if __name__ == '__main__':
    # smi = 'COC1=CC=CC([C@](O2)(CN3C=CN=C3)OC[C@H]2COC4=CC=CC=C4)=C1'
    smi = 'COC1CC(N)C1'
    # smi = 'CO'

    xyz = smi2xyz(smi)
    # print(xyz)

    shutil.rmtree(TMP_DIR, ignore_errors=True)
    os.makedirs(TMP_DIR, exist_ok=True)
    os.chdir(TMP_DIR)

    xyz_fname = 'struc_01.xyz'
    with open(xyz_fname, 'w') as f:
        f.write(xyz)

    xtb_optimize(xyz_fname, method="gff")

    # if os.path.exists('xtbopt.xyz'):
    #     print('Success geometry optimization')

    xyz_fname = 'xtbopt.xyz'
    conformer_gen(
        xyz_fname,
        method="gfnff",
        ewin=8.0,
        mdlen=5.0,
        n_jobs=8,
    )

    xyz_fname = 'crest_conformers.xyz'
    conformer_screen(
        xyz_fname,
        method="gfn2",
        ewin=12.0,
        n_jobs=16,
        verbose=False)


    if os.path.exists('crest_ensemble.xyz'):
        print('Success conformer screen')

    shutil.copy(f'crest_best.xyz', '../../data')