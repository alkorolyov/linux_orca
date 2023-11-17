import logging

logging.basicConfig(level=logging.DEBUG)

if __name__ == '__main__':
    import os
    import molli as ml

    from indigo import Indigo
    indigo = Indigo()

    from openbabel import openbabel
    obmol = openbabel.OBMol()
    obconv = openbabel.OBConversion()

    def smi_to_mol2(smiles):
        mol = indigo.loadMolecule(smiles)
        mol.layout()
        obconv.SetInAndOutFormats("mol", "mol2")
        obconv.ReadString(obmol, mol.molfile())
        return obconv.WriteString(obmol)

    # smi = 'CC(=O)N1CCN(CC1)C2=CC=C(C=C2)OC[C@@H]3CO[C@@](O3)(CN4C=CN=C4)C5=C(C=C(C=C5)Cl)Cl'
    smi = 'BrC1=CC=C(N(CCC)CCC)C=C1'

    molli_mol = ml.Molecule.from_mol2(smi_to_mol2(smi))
    col = ml.Collection(name="molecule", molecules=[molli_mol])

    tmp_dir = "tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    crest = ml.CRESTDriver(
        name="confs",
        scratch_dir=tmp_dir + "/crest_scratch_1/",
        nprocs=8,
    )

    concur_1 = ml.Concurrent(
        col,
        backup_dir=tmp_dir + "/crest_search/",
        logfile=tmp_dir + "/out1.log",
        update=60,
        # timeout=10000,
        timeout=600,
        concurrent=1,
    )
    
    print("conformer search beginning")
    output = concur_1(crest.conformer_search)(ewin=8, mdlen=5, constr_val_angles=[])
    print("searched conf\n", output)
    buffer = []
    tracking = {}  # Used to track progress
    for i, k in enumerate(output):
        if isinstance(k, ml.Molecule):
            tracking[k.name] = True
            buffer.append(k)
        else:
            tracking[col.molecules[i].name] = False
    print(buffer)
    if len(buffer) == 0:
        raise Exception(f"Error calculating conformers - search step failed")

    col2 = ml.Collection(
        name="searched", molecules=buffer
    )  # These have undergone a conformer search.
    print(col2.molecules)
    concur_2 = ml.Concurrent(
        col2,
        backup_dir=tmp_dir + "/crest_screen/",
        logfile=tmp_dir + "/out2.log",
        update=30,
        # timeout=10000,
        timeout=600,
        concurrent=1,
    )
    crest = ml.CRESTDriver(
        name="confs",
        scratch_dir=tmp_dir + "/crest_scratch_2/",
        nprocs=8,
    )
    print("conformer screen beginning")
    output2 = concur_2(crest.confomer_screen)(
        method="gfn2", ewin=12
    )  # These get screened to prune out unreasonable structures and reopt.
    buffer2 = []
    print("screened conf\n", output2)

    assert len(output2) > 0
    for j, k in enumerate(output2):
        if isinstance(k, ml.Molecule):  # By definition, must have succeeded before
            buffer2.append(k)
        else:  # Failed - make sure it is set to False in tracking dictionary
            if (
                    tracking[col2.molecules[j].name] == True
            ):  # If it worked before, but failed here, set it to False.
                tracking[col2.molecules[j].name] = False
            else:
                pass  # Already set to False, can skip.
    if len(buffer2) == 0:
        raise Exception(f"Error calculating conformers - screen step failed")
    col3 = ml.Collection(name="screened", molecules=buffer2)
    assert len(buffer2) > 0
    conformers = col3
    failures = []
    for key, val in tracking.items():
        if val == False:
            failures.append(col[key])
        elif val == True:
            pass
    if (
            len(failures) > 0
    ):  # If molecules failed, warn the user. Put this in a log or something later. DEV
        raise Warning(
            f"Molecules failed conformer search: {[f.name for f in failures]}"
        )
    failures = failures
    if len(failures) > 0:
        towrite = ml.Collection(name="failed", molecules=failures)
        towrite.to_zip(
            tmp_dir + "/input_mols_failed_conf_step.zip"
        )
    assert len(conformers.molecules) > 0
    print("Script finished")

