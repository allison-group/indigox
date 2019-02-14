#!/usr/bin/env python3

import indigox as ix
from pathlib import Path

# You always need a forcefield
ff = ix.GenerateGROMOS54A7()

# Generate and save the Athenaeums used in the paper
# On our system this takes approximately 6.5 minutes to execute and uses a maximum of 3.5GB of memory.
# The saved files are 463 kB and 684.2 MB respectively.
def LoadAndSaveAthenaeums():
    settings = ix.Athenaeum.Settings
    man_ath = ix.Athenaeum(ff)
    man_ath.SetBool(settings.SelfConsistent)
    auto_ath = ix.Athenaeum(ff, 1)

    # Need to set a larger than default limit as some of the amino acid molecules are large
    auto_ath.SetInt(settings.MoleculeSizeLimit, 60)
    
    # Load all the molecules
    for aa in Path("SourceMolecules").glob("*.frag"):
        mol, fragments = ix.LoadFragmentFile(aa, ff)
        mol.SetName(aa.stem)
        # add the fragments to the manual Athenaeum
        for frag in fragments: man_ath.AddFragment(frag)
        # add the molecule to the automatic Athenaeum
        auto_ath.AddAllFragments(mol)

    # Save the athenaeums
    ix.SaveAthenaeum(man_ath, "ManualAthenaeum.ath")
    ix.SaveAthenaeum(auto_ath, "AutomaticAthenaeum.ath")

# Load athenaeums and run the cherrypicker algorithm on the 3 molecules used in the paper
# On our system, loading the Athenaeums takes about 12 seconds, then running Loading, Running CherryPicker, and Saving the test molecules takes about 5 seconds.
def RunCherryPicker():
    settings = ix.algorithm.CherryPicker.Settings
    # Load the athenaeums
    man_ath = ix.LoadAthenaeum("ManualAthenaeum.ath")
    auto_ath = ix.LoadAthenaeum("AutomaticAthenaeum.ath")  # Is a large file so may take a few seconds to load

    # Create the CherryPicker algorithm to use
    cherrypicker = ix.algorithm.CherryPicker(ff)
    cherrypicker.AddAthenaeum(man_ath)  # add the manual one first so it's used first
    cherrypicker.AddAthenaeum(auto_ath)
    
    # Set the CherryPicker options we want
    cherrypicker.SetInt(settings.MinimumFragmentSize, 2)
    cherrypicker.SetInt(settings.MaximumFragmentSize, 20)

    # Load each of the test molecules and run cherrypicker
    for test in Path("TestMolecules").glob("*.pdb"):
        mol = ix.LoadPDBFile(test, test.with_suffix(".ixd"))
        parameterised = cherrypicker.ParameteriseMolecule(mol)

        # save the parameterisation in ITP format
        ix.SaveITPFile(test.with_suffix(".itp"), mol, parameterised)

if __name__ == "__main__":
    LoadAndSaveAthenaeums()
    RunCherryPicker()
