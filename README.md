# indigox

Package containing the CherryPicker algorithm. It is intended that indigox is used primarily as a Python module.

## Requirements
- CMake version >= 3.7
- C++17 compilent compiler. In our tests this means GCC >= 7
- [optional] Boost. The required libraries are provided from boost version 1.66.0
- [optional] Python >= 3.4, with Python dev tools installed, for building the python bindings. Highly recommended.
- [optional] doxygen for building the documentation

## Installation
Tested on macOS and Linux. If you're on Windows, you're on your own. We have not tested on a Windows system.

To obtain everything for compilation, you will need to clone recursively, that is:

git clone --recurse-submodules https://github.com/allison-group/indigox.git

mkdir build && cd build

cmake .. && make

make install  (depending on where your python packages are installed, this may need sudo)

Examples are provided in the examples directory, with a worked explanation of the code in the [tutorial](https://allison-group.github.io/indigox/tutorials.html). Documentation can be built by make doc, or a full version is available online at https://allison-group.github.io/indigox/

If you are developing using the python examples, note that a full build and reinstallation of indigox is needed for the example to pick up on changes. The example relies on the installed indigox version. 

## Python Bindings

The easiest way to use Cherry Picker is through the indigox python module. This will be built and installed on your system by following the instructions above.

Some notes on the bindings:
 * Installation requires the python dev tools. 
 ```
sudo apt-get install python3-dev 

or

sudo yum install python3-devel
```
 * It will be installed in the site-packages location of your default python executable. You can install in specific environments by changing the PYTHON_SITE_PACKAGES variable in [`indigox/src/python/CMakeLists.txt`](./src/python/CMakeLists.txt). 
 * The python library has most of the file format input and output methods. See [parser.py](./src/python/io/parser.py) and [serialiser.py](./src/python/io/serialiser.py) for some useful methods. 

## File formats

Indigox uses two new file formats, and supports a series of other ones. 

### IXD format

This format stores formal charges, implicit hydrogen counts, bond orders, and the total molecular charge of a molecule. It has three entry types which must be on their own line but can appear in any order:

- ATOM entries. ```ATOM  X  FC  H``` 'ATOM' is the line header, and must always be included. X is the atom index, FC is the formal charge, and H is the number of implicit hydrogens.
- BOND entries. ```BOND  X  Y  O``` 'BOND' is the line header, and must always be included. X and Y are the atom indices, and O is the order. O can be 1 - 4 for single to quadruple bonds, 5 for aromatic, 6 for 1.5 and 7 for 2.5.
- MOLECULE entry. ```MOLECULE  Q``` 'MOLECULE' is the line header, and must always be included. Q is the total molecular charge. This must match the sum of formal charges otherwise an error will be thrown.

Note that not every atom/bond needs to be specified in an IXD file, but each line entry expects every value to be included for that atom/bond. 

### Fragment files

This format is used to describe fragments of larger molecules to include when making an Athenaeum. By manually specifying the fragments, you have greater control over which pieces of molecules are recognised and parameterised as a unit. Our examples, for instance, include fragment definitions for full amino acids, allowing them to be matched and parameterised as a whole. 

The fragment file has three types of information blocks:
- The 'MOLECULE' block. Must always come first in the file, and must only occur once. It specifies a coordinate file (.pdb), parameter file (.mtb or .itp), and IXD file, from which these fragments are built. These paths must all be included, must be on separate lines, and are specified relative to the fragment file.
- The 'FRAGMENT' block. Contains the list of atom indices to use as a core region for a fragment. 
- The 'OVERLAP' block. Contains the list of atom indices to use as the overlap region for the preceding fragment.

A simple example (more can be found in examples->CherryPicker->SourceMolecules):

```
MOLECULE
Alanine_ua.pdb
Alanine.top
Alanine.ixd
END
FRAGMENT
24 25 26 27 28 29
END
OVERLAP
22 23 30 31
END
FRAGMENT
1 2 3 4 5 6 7 8 9 10
11 12 13 14 15 16 17 18 19 20 
21 22 23
END
OVERLAP
END
```

A fragment file can have any number of FRAGMENT blocks, but each must be followed by an OVERLAP block, even if the OVERLAP block has no atoms. 

### Additional file formats

CherryPicker can read the following file formats:
- PDB - to load a molecule. CherryPicker expects PDB with CONECT records. 
- ITP and MTB - to load parameters for a molecule. Usually read in next to the PDB. ITP parsing assumes a GROMOS-style forcefield.
- IXD - bespoke format that can be read in next to ITP or MTB to provide further parameters to a molecule.
- IFP - to read in forcefield parameters. Currently only supports GROMOS style IFP files, but CherryPicker comes with a pre-defined implementation of GROMOS 54A7.

And it can output the following formats:
- PDB - for molecular structure.
- ITP, RTP, IXD - for outputting parameters of the molecule.

## Options

Indigox exposes a number of customisation and optimisation options. To set options, see the example python below. Note that the Athenaeum and CherryPicker settings are separate.

```
import indigox as ix

force_field = ix.GenerateGROMOS54A7() # always need a force field first

# Settings for Athenaeum generation
ath_settings = ix.Athenaeum.Settings
example_ath = ix.Athenaeum(force_field)

example_ath.SetInt(ath_settings.MoleculeSizeLimit, 60)
example_ath.SetBool(ath_settings.SelfConsistent)

# Settings for Cherry Picker
cp_settings = ix.algorithm.CherryPicker.Settings
cherrypicker = ix.algorithm.CherryPicker(force_field)

# Set fragment size boundaries
cherrypicker.SetInt(settings.MinimumFragmentSize, 2)
cherrypicker.SetInt(settings.MaximumFragmentSize, 20)
```
For available Athenaeum generation options, see the [settings section of the Athenaeum class definition](https://allison-group.github.io/indigox/classindigox_1_1_athenaeum.html#a4c05250977b08f4d2fe1fe19dd1df0fa). 
For available CherryPicker options, see the [settings section of the CherryPicker class definition](https://allison-group.github.io/indigox/classindigox_1_1algorithm_1_1_cherry_picker.html#a4c05250977b08f4d2fe1fe19dd1df0fa).

### Formal Charge and Bond Order assignment

If the CherryPicker setting ```CalculateElectrons``` is true (it's false by default), CherryPicker will run the indigo-bondorder electron assignment algorithm on your molecule first. The algorithm is detailed [here](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0340-0) and calculates bond orders and formal charges.

Use the CherryPicker option ```ElectronMethod``` to set the underlying search algorithm:
- 0 = Local Optimisation
- 1 = A*
- 2 = FPT (Fixed Parameter Tractable) (the default setting, as the most accurate method)

If multiple resonance structures are found the program will pause and wait for user input to choose between them. To get around this, set ```NoInput``` to true. The first resonance structure found will be used. NOTE: We can't guarantee anything about this struture, so this setting is only recommended for very simple molecules that need to be processed in batches.

This version of indigo-bondorder uses a set of default options (outlined below). If you need other settings, use the standalone algorithm located [here](https://github.com/allison-group/indigo-bondorder) to parameterise and export your molecule before importing it into CherryPicker.

##### Default Formal Charge and Bond Order options

Option name | Value | Note
------ | :------: | ------------
USE_ELECTRON_PAIRS | true | Calculate electrons as pairs to reduce search space
PREPLACE_ELECTRONS | true | Preplace 6 electrons on halogens singly bonded to carbons to reduce search space
ALLOWED_ELEMENTS | H, C, N, O, S, P, F, Cl, Br | Encountering elements outside this list will cause electron assignment to fail
MAXIMUM_BOND_ORDER | 3 | 