#!/usr/bin/env python3

import indigox as ix
from pathlib import Path
import time
import fnmatch
import os

t_ff = time.time()

# You always need a forcefield
ff = ix.GenerateGROMOS54A7()
total_ff = time.time() - t_ff

manualAthPath = "ManualAthenaeum.ath"
autoAthPath = "AutomaticAthenaeum.ath"

def LoadAndSaveAthenaeums():
    """ Generate and save the Athenaeums used in the paper
    On our system this takes approximately 6.5 minutes to execute and uses a maximum of 3.5GB of memory.
    The saved files are 463 kB and 684.2 MB respectively.
    :return:
    """
    pre_init_aths = time.time_ns()

    settings = ix.Athenaeum.Settings
    man_ath = ix.Athenaeum(ff)
    man_ath.SetBool(settings.SelfConsistent)
    auto_ath = ix.Athenaeum(ff, 1)

    # Need to set a larger than default limit as some of the amino acid molecules are large
    auto_ath.SetInt(settings.MoleculeSizeLimit, 60)

    print("Initialised Athenaeums in %d ns" % (time.time_ns() - pre_init_aths))

    mol_path = "SourceMolecules"
    mol_extn = "*.frag"

    index = 1
    mol_count = len(fnmatch.filter(os.listdir(mol_path), mol_extn))
    load_times, man_ath_times, auto_ath_times = [], [], []

    print("Loading molecules into Athenaeums...")
    pre_loading = time.time_ns()

    # Load all the molecules
    for aa in Path(mol_path).glob(mol_extn):
        print("\tLoading Amino acid #{} of {}".format(index, mol_count))
        index += 1

        pre_loadfrag = time.time_ns()

        mol, fragments = ix.LoadFragmentFile(aa, ff)
        mol.SetName(aa.stem)

        post_loadfrag = time.time_ns()

        # add the fragments to the manual Athenaeum
        for frag in fragments:
            man_ath.AddFragment(frag)
        post_man_ath = time.time_ns()

        # add the molecule to the automatic Athenaeum
        auto_ath.AddAllFragments(mol)
        post_auto_ath = time.time_ns()

        load_times.append(post_loadfrag - pre_loadfrag)
        man_ath_times.append(post_man_ath - post_loadfrag)
        auto_ath_times.append(post_auto_ath - post_man_ath)

    print("\nFinished loading Athenaeums in {:10.3f} s".format((time.time_ns() - pre_loading)/(10**9)))
    print("\tTotal loading from file: {:10.3f} s".format(sum(load_times)/(10**9)))
    print("\tTotal Manual Athenaeum:\t {:10.3f} s".format(sum(man_ath_times)/(10**9)))
    print("\tTotal Automatic Athenaeum: {:10.3f} s".format(sum(auto_ath_times)/(10**9)))

    pre_save = time.time_ns()
    # Save the athenaeums
    ix.SaveAthenaeum(man_ath, manualAthPath)
    ix.SaveAthenaeum(auto_ath, autoAthPath)
    print("Saved Athenaeums to file in  {:10.3f} s\n".format((time.time_ns() - pre_save)/(10**9)))

def RunCherryPicker():
    """ Load athenaeums and run the cherrypicker algorithm on the 3 molecules used in the paper
    On our system, loading the Athenaeums takes about 12 seconds,
    then running Loading, Running CherryPicker, and Saving the test molecules takes about 5 seconds.
    :return:
    """

    settings = ix.algorithm.CherryPicker.Settings
    # Load the athenaeums
    man_ath = ix.LoadAthenaeum(manualAthPath)
    auto_ath = ix.LoadAthenaeum(autoAthPath)  # Is a large file so may take a few seconds to load

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
    before_loadAth = time.time()
    if not(Path(manualAthPath).is_file() and Path(autoAthPath).is_file()):
        LoadAndSaveAthenaeums()

    before_CP = time.time()
    RunCherryPicker()

    print("\nForce field: {:10.3f} s".format(total_ff))
    print("Load Athenaeums: {:10.3f} s".format(before_CP - before_loadAth))
    print("Run CherryPicker: {:10.3f} s".format(time.time() - before_CP))
