#!/usr/bin/env python3

"""
Construct Athenaeum and Parameterise molecule example.

This example creates two Athenaeums out of the amino acid files in ./SourceMolecules:
- AutomaticAthenaem: contains all possible fragments of the amino acids, with 1 atom
overlaps
- ManualAthenaeum: contains fragments as defined by the .frag files in ./SourceMolecules.

This example runs CherryPicker using these Athenaeums to parameterise the molecules
in ./TestMolecules. The parameters are saved in the same folder in .itp format.

Note that ManualAthenaeum is added to CherryPicker first, meaning full amino acids will
be preferentially identified before the algorithm resorts to smaller fragments in the
AutomaticAthenaeum.
"""

import indigox as ix
from pathlib import Path
import time
import fnmatch
import os

# You always need a forcefield
ff = ix.GenerateGROMOS54A7()

manualAthPath = "ManualAthenaeum.ath"
autoAthPath = "AutomaticAthenaeum.ath"


def LoadAndSaveAthenaeums():
    """ Generate and save the Athenaeums used in the paper
    On our system this takes approximately 6.5 minutes to execute and uses a maximum of 3.5GB of memory.
    The saved files are 463 kB and 684.2 MB respectively.
    :return:
    """

    settings = ix.Athenaeum.Settings
    man_ath = ix.Athenaeum(ff)
    man_ath.SetBool(settings.SelfConsistent)
    auto_ath = ix.Athenaeum(ff, 1)

    # Need to set a larger than default limit as some of the amino acid molecules are large
    auto_ath.SetInt(settings.MoleculeSizeLimit, 60)

    mol_path = "SourceMolecules"
    mol_extn = "*.frag"

    index = 1
    mol_count = len(fnmatch.filter(os.listdir(mol_path), mol_extn))

    print("Loading molecules into Athenaeums...")

    # Load all the molecules
    for aa in Path(mol_path).glob(mol_extn):
        print("\tLoading Amino acid #{} of {}".format(index, mol_count))
        index += 1

        mol, fragments = ix.LoadFragmentFile(aa, ff)
        mol.SetName(aa.stem)

        # add the fragments to the manual Athenaeum
        for frag in fragments:
            man_ath.AddFragment(frag)

        # add the molecule to the automatic Athenaeum
        auto_ath.AddAllFragments(mol)

    # Save the athenaeums
    ix.SaveAthenaeum(man_ath, manualAthPath)
    ix.SaveAthenaeum(auto_ath, autoAthPath)


def RunCherryPicker():
    """ Load athenaeums and run the cherrypicker algorithm on the 3 molecules used in the paper.
    On our system, loading the Athenaeums takes about 12 seconds,
    then running Loading, Running CherryPicker, and Saving the test molecules takes about 5 seconds.
    :return:
    """

    # Load the athenaeums
    man_ath = ix.LoadAthenaeum(manualAthPath)
    auto_ath = ix.LoadAthenaeum(autoAthPath)  # Is a large file so may take a few seconds to load

    # Create the CherryPicker algorithm to use
    cherrypicker = ix.algorithm.CherryPicker(ff)
    cherrypicker.AddAthenaeum(man_ath)  # add the manual one first so it's used first
    cherrypicker.AddAthenaeum(auto_ath)

    # Set the CherryPicker options we want
    settings = ix.algorithm.CherryPicker.Settings
    cherrypicker.SetInt(settings.MinimumFragmentSize, 2)
    cherrypicker.SetInt(settings.MaximumFragmentSize, 20)

    # Load each of the test molecules and run cherrypicker
    for test in Path("TestMolecules").glob("*.pdb"):
        mol = ix.LoadPDBFile(test, test.with_suffix(".ixd"))
        parameterised = cherrypicker.ParameteriseMolecule(mol)

        # save the parameterisation in ITP format
        ix.SaveITPFile(test.with_suffix(".itp"), mol, parameterised)
        print()


if __name__ == "__main__":
    before_loadAth = time.time()
    if not (Path(manualAthPath).is_file() and Path(autoAthPath).is_file()):
        LoadAndSaveAthenaeums()

    before_CP = time.time()
    RunCherryPicker()

    print("\nLoad Athenaeums: {:10.3f} s".format(before_CP - before_loadAth))
    print("Run CherryPicker: {:10.3f} s".format(time.time() - before_CP))
