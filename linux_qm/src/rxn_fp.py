import logging
import numpy as np
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdMolDescriptors


def fp2numpy(x):
    return np.frombuffer(x.ToBitString().encode(), 'u1') - ord('0')


def get_morgan_fp(input, radius: int = 2, bits: int=1024, **kwargs):
    try:
        mol = None
        if isinstance(input, str):
            mol = Chem.MolFromSmiles(input)
        elif isinstance(input, Chem.rdchem.Mol):
            mol = input
        else:
            raise ValueError('Fingerprint Error: Expected rdmol or smiles string')
        # print(kwargs)
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=bits, **kwargs)
        return fp2numpy(fp)
    except Exception as e:
        logging.debug(f"{type(e).__name__} during generating fingerprint: {e}")


def init_mol(mol):
    mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(mol)


@dataclass
class Bond:
    begin: int
    end: int
    order: Chem.rdchem.BondType


@dataclass
class Reaction:
    ephile: str
    nphile: str
    product: str


def get_reacting_bonds(rxn, reacting_atoms):
    print(reacting_atoms)
    reacting_bonds = []
    for ridx, raids in enumerate(reacting_atoms):
        reacting_bonds.append([])
        reactant = rxn.GetReactantTemplate(ridx)
        N = len(raids)
        for i in range(N - 1):
            for j in range(i + 1, N):
                b = reactant.GetBondBetweenAtoms(raids[i], raids[j])
                if b is None:
                    continue
                bond = Bond(raids[i], raids[j], b.GetBondType())
                reacting_bonds[ridx].append(bond)
    return reacting_bonds


def get_bonds(atoms_ids, mol):
    bonds = []
    N = len(atoms_ids)
    for i in range(N - 1):
        for j in range(i + 1, N):
            b = mol.GetBondBetweenAtoms(atoms_ids[i], atoms_ids[j])
            if b is None:
                continue
            bond = Bond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondType())
            bonds.append(bond)
    return bonds


def get_substruct(atoms_ids, mol):
    sub = Chem.RWMol()
    # mapping dictionary
    # when atoms added to the new mol, they new indices
    # old_atom_idx -> new_atom_idx
    atom_map = {}
    for aidx in atoms_ids:
        a_num = mol.GetAtomWithIdx(aidx).GetAtomicNum()
        # print(mol.GetAtomWithIdx(aidx).GetPropsAsDict())
        atom_map[aidx] = sub.AddAtom(Chem.Atom(a_num))
    for bond in get_bonds(atoms_ids, mol):
        sub.AddBond(atom_map[bond.begin], atom_map[bond.end], bond.order)

    mol = sub.GetMol()
    mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(mol)

    return mol


def get_neighbours(atoms_ids, mol):
    nei_ids = []
    for aidx in atoms_ids:
        a = mol.GetAtomWithIdx(aidx)
        for nei in a.GetNeighbors():
            nei_idx = nei.GetIdx()
            # print('nei_idx', nei_idx)
            if nei_idx not in nei_ids and nei_idx not in atoms_ids:
                nei_ids.append(nei_idx)
    return nei_ids


def get_layers(orig_ids, mol, max_dist=2):
    """
    Generate a list of incremental layers. Each one contains an incrementing
    list of atom indices each with increasing distance from original set.
    The step is 1 bond distance until max_radius bond distance reached.
    :param atoms_ids: list of starting atoms indices ex. [0]
    :param mol: parent mol
    :param max_dist: max topological radius in bonds
    :return: list of increasing layers, ex. '[[0], [0, 1, 4], [0, 1, 4, 2, 3, 5]]'
    """
    r = 0
    atoms_ids = list(orig_ids)
    res = [atoms_ids]
    while r < max_dist:
        tmp = res[-1]
        res.append(tmp.copy() + get_neighbours(tmp, mol))
        r += 1
    return res


def get_layers_mols(atoms_ids, mol, max_dist=2):
    layers_aids = get_layers(atoms_ids, mol, max_dist=max_dist)
    return [get_substruct(aids, mol) for aids in layers_aids]


def map_react_atoms_to_product(reacting_atoms, rxn):
    """
    Map reacting atoms idx to the corresponding idx in product
    :param reacting_atoms: tuple of reacting atom idx
    :param rxn: rdkit Reaction
    :return: list of reacting atom idx in product
    """

    map_nums = []
    for ridx, atoms_ids in enumerate(reacting_atoms):
        reactant = rxn.GetReactantTemplate(ridx)
        for aidx in atoms_ids:
            a = reactant.GetAtomWithIdx(aidx)
            map_num = a.GetAtomMapNum()
            if map_num:
                map_nums.append(map_num)
    res = []
    for a in rxn.GetProducts()[0].GetAtoms():
        if a.GetAtomMapNum() in map_nums:
            res.append(a.GetIdx())
    return res


class RxnFpGenerator:
    """
    Multi layer Reaction fingerprint:
    Three layers:
     - full fp for components:
            ephile                          256 bit
            nphile                          256 bit
            product                         512 bit
     - layers for reacting center:
        dist 0 (only reacting atoms)
            ephile                          128  bit
            nphile                          128  bit
            product                         256 bit
        dist 1
            ephile                          128  bit
            nphile                          128  bit
            product                         256 bit
        dist 2
            ephile                          128  bit
            nphile                          128  bit
            product                         256 bit

        Total                               2048 bit
    """

    def __init__(self,
                 layers_shape=(
                     [256, 256, 512],    # full fp for rxn reagents and product
                     [  0,   0,   0],      # dist 0 layer (only reacting atoms)
                     [128, 128, 256],    # dist 1 layer
                     [128, 128, 256],    # dist 2 layer
                 ),
                 fp_radius=2,
                 calc_reags=True,
                 calc_layers=True):
        self.calc_reags = calc_reags
        self.calc_layers = calc_layers
        self.fp_radius = fp_radius
        self.max_layers_dist = len(layers_shape) - 1
        self.layers_shape = layers_shape

    def _calc_init(self, rxn_smi):
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
        rxn.Initialize()
        ephile, nphile = rxn.GetReactants()
        rxn_prod = rxn.GetProducts()[0]
        init_mol(ephile)
        init_mol(nphile)
        init_mol(rxn_prod)

        self.reacting_atoms = rxn.GetReactingAtoms()
        self.rxn = rxn
        self.ephile = ephile
        self.nphile = nphile
        self.rxn_prod = rxn_prod

        self.reag_fp = np.array([])
        self.layers_fp = np.array([])

    def _calc_reag_fp(self):
        """
        Calculates stacked fp for full reagents - ephile, nphile, rxn_prod
        @param layer_size: array of bit sizes, ex. '[128, 128, 256]'
        @return: numpy array of stacked fp's
        """
        e_fp = get_morgan_fp(self.ephile, self.fp_radius, self.layers_shape[0][0])
        n_fp = get_morgan_fp(self.nphile, self.fp_radius, self.layers_shape[0][1])
        p_fp = get_morgan_fp(self.rxn_prod, self.fp_radius, self.layers_shape[0][2])
        fp_list = [np.empty(0, dtype='uint8') if fp is None else fp for fp in [e_fp, n_fp, p_fp]]
        self.reag_fp = np.hstack(fp_list)

    def _calc_layers(self):
        """
        Calculate layers per each reagent and product.
        @return: (dict) python dictionary of layers
        """
        reacting_atoms_prod = map_react_atoms_to_product(self.reacting_atoms, self.rxn)
        max_dist = self.max_layers_dist

        layers_dict = {}

        layers_dict['ephile'] = get_layers(self.reacting_atoms[0], self.ephile, max_dist=max_dist)
        layers_dict['nphile'] = get_layers(self.reacting_atoms[1], self.nphile, max_dist=max_dist)
        layers_dict['rxn_prod'] = get_layers(reacting_atoms_prod, self.rxn_prod, max_dist=max_dist)
        self.layers_dict = layers_dict

    def _calc_layers_fp(self):
        """
        Calculate and stack fp for each layer.
        @return: stacked fp as numpy array
        """
        self._calc_layers()

        fps = []
        # print(self.layers_shape)
        # print(self.max_layers_dist)
        for i in range(1, len(self.layers_shape)):
            fp_e = get_morgan_fp(self.ephile, self.fp_radius, self.layers_shape[i][0], fromAtoms=self.layers_dict['ephile'][i])
            fp_n = get_morgan_fp(self.nphile, self.fp_radius, self.layers_shape[i][1], fromAtoms=self.layers_dict['nphile'][i])
            fp_p = get_morgan_fp(self.rxn_prod, self.fp_radius, self.layers_shape[i][2], fromAtoms=self.layers_dict['rxn_prod'][i])
            fp_list = [np.empty(0, dtype='uint8') if fp is None else fp for fp in [fp_e, fp_n, fp_p]]
            fps.append(np.hstack(fp_list))
            # print('fps.shape', np.hstack([fp_e, fp_n, fp_p]).shape)

        # print('len(fps)', len(fps))
        # print('fps[0]', fps[0])
        self.layers_fp = np.hstack(fps)

    def get_rxnfp(self, rxn_smi):
        self._calc_init(rxn_smi)
        if self.calc_reags:
            self._calc_reag_fp()
        if self.calc_layers and len(self.layewrs_shape) > 1:
            self._calc_layers_fp()
        return np.hstack([self.reag_fp, self.layers_fp])


