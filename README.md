# indigox

Package containing the CherryPicker algorithm. It is intended that indigox is used primarily as a Python module.

## Requirements
- CMake version >= 3.7
- C++17 compilent compiler. In our tests this means GCC >= 7
- [optional] Boost. The required libraries are provided from boost version 1.66.0
- [optional] Python >= 3.4, for building the python bindings. Highly recommended
- [optional] doxygen for building the documentation

## Installation
Tested on macOS and Linux. If you're on Windows, you're on your own. We have no means of testing on a Windows system.

To correctly obtain everything for compilation, you will need to clone recursively, that is:

git clone --recurse-submodules https://github.com/allison-group/indigox.git

mkdir build && cd build

cmake .. && make

make install  (depending on where your python packages are installed, this may need sudo)

Examples are provided in the examples directory, with a worked explanation of the code in the [tutorial](https://allison-group.github.io/indigox/tutorials.html). Documentation can be built by make doc, or a full version is available online at https://allison-group.github.io/indigox/

If you are developing using the python examples, note that a full build and reinstallation of indigox is needed for the example to pick up on changes. The example relies on the installed indigox version. 

## File formats

We use two new file formats. The first is the IXD format. This format allows for formal charges, implicit hydrogen, bond orders and total molecular charges to be defined. It has three possible entry types:
- "ATOM  X  FC  H". This type for specifing things to do with atoms. The header 'ATOM' must always be included. X is the atom index which the settings should be applied to, FC is the formal charge and H is the number of implicit hydrogens.
- "BOND  X  Y  O". This type is for specifing bond order between two atoms. The header 'BOND' must always be included. X and Y are the atom indices and O is the order. O can be 1 - 4 for single to quadruple bonds, 5 for aromatic, 6 for 1.5 and 7 for 2.5.
- "MOLECULE  Q". This type is for specifing the total molecular charge of a molecule. The header 'MOLECULE' must always be included. Q is the total molecular charge. This must match the sum of formal charges otherwise an error will be thrown.

The second file format is the fragment definition file. This format contains three blocks, each opened with a header and terminated with 'END'.
- The 'MOLECULE' block must always come first in the file, and must only occur once. It contains three lines indicating which files to load for the molecule definition. First, a coordinate file, second a parameter file and last an IXD file. These file paths are relative to the fragment file path.
- The 'FRAGMENT' block contains the list of atom indices to use for the core region in fragment generation. This is a free-form block taking any number of lines with any number of indices per line.
- The 'OVERLAP' block contains the list of atom inidices to use for the overlap region in fragment generation. Like the 'FRAGMENT' block, this is free form. Each occurance of a 'FRAGMENT' block requires the next block to be an 'OVERLAP' block, even if the 'OVERLAP' block is empty.

## Options

Indigox exposes a number of customisation and optimisation options. To set options using the python module:

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

### Available options

#### Athenaeum generation:


#### Molecule parameterisation


#### Formal Charge and Bond Order assignment

If CalculateElectrons is set to true, CherryPicker will run the indigo-bondorder electron assignment algorithm detailed [here](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0340-0) to calculate bond orders and formal charges, before running the CherryPicker algorithm.

Use the option ElectronMethod to set the underlying search algorithm:
0 = Local Optimisation
1 = A*
2 = FPT (Fixed Parameter Tractable)

By default FPT is used and recommended as the most accurate method. If time is the limiting resource and accuracy can be sacrificed, Local Optimisation is the fastest.

This implementation uses a set of default options outlined below. If you need to use other settings, you can use the standalone algorithm to parameterise and export your molecule first. It has a more accessible python binding library for its options. Find the source code [here](https://github.com/allison-group/indigo-bondorder).

##### Default Formal Charge and Bond Order options

Option name | Value | Note
------ | :------: | ------------
USE_ELECTRON_PAIRS | true | Calculate electrons as pairs to reduce search space
PREPLACE_ELECTRONS | true | Preplace 6 electrons on halogens singly bonded to carbons to reduce search space
ALLOWED_ELEMENTS | H, C, N, O, S, P, F, Cl, Br | Encountering elements outside this list will cause electron assignment to fail
MAXIMUM_BOND_ORDER | 3 | 