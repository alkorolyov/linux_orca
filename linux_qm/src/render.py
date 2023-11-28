
from rdkit.Chem import rdChemReactions
from indigo import Indigo
from indigo.renderer import IndigoRenderer

indigo = Indigo()
renderer = IndigoRenderer(indigo)

indigo.setOption("render-output-format", "png")
# indigo.setOption("render-superatom-mode", "collapse")
indigo.setOption("render-coloring", True)
# indigo.setOption("render-base-color", "1, 1, 1")
indigo.setOption("render-relative-thickness", 1.5)
indigo.setOption("render-highlight-thickness-enabled", True)
indigo.setOption("render-highlighted-labels-visible", True)
# indigo.setOption("render-highlight-color", "1, 0.4, 0")

from IPython.display import display, Image
from indigo.renderer import IndigoRenderer

def draw_indigo_obj(obj, size=None):
    if size:
        if isinstance(size, int):
            size = f"{size}, {size}"
        indigo.setOption("render-image-size", size)
    display(Image(renderer.renderToBuffer(obj)))


def draw_reaction_smi(rxn_smi: str, highlight: list = None, product_only=False):
    # auto size
    indigo.setOption("render-image-size", "-1, -1")

    try:
        rxn = indigo.loadReaction(rxn_smi)

        if highlight:
            indigo.setOption("render-coloring", False)
            for mol in rxn.iterateMolecules():
                for atom in mol.iterateAtoms():
                    map_num = rxn.atomMappingNumber(atom)
                    if map_num in highlight:
                        atom.highlight()
        if product_only:
            indigo.setOption("render-image-size", "300, 300")
            for prod in rxn.iterateProducts():
                draw_indigo_obj(prod, )
        else:
            draw_indigo_obj(rxn)
        indigo.setOption("render-coloring", True)
    except Exception as e:
        print(f"Parsing error: {e}\nSMILES: {rxn_smi}")


def draw_reacting_mapnums(rxn_smi, verbose=True):
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
    rxn.Initialize()

    # reacting_atoms = rxn.GetReactingAtoms(mappedAtomsOnly=True)
    reacting_atoms = rxn.GetReactingAtoms()
    if verbose:
        print('Reacting Atom Idx:', reacting_atoms)
    reacting_mapnums = []

    for ridx, reacting in enumerate(reacting_atoms):
        reactant = rxn.GetReactantTemplate(ridx)
        for raidx in reacting:
            atom = reactant.GetAtomWithIdx(raidx)
            mapnum = atom.GetAtomMapNum()
            reacting_mapnums.append(mapnum)
            # print(Chem.MolToSmiles(reactant))
            if verbose:
                print('Mapped reacting atom:', atom.GetSymbol(), mapnum)
    if verbose:
        draw_reaction_smi(rxn_smi, highlight=reacting_mapnums)
    return reacting_mapnums