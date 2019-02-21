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

make install

Examples are provided in the examples directory. Documentation can be built by make doc. Alternatively, a version of the documentation up to date with the most recent release is available at https://allison-group.github.io/indigox/

## File formats

We use two new file formats. The first is the IXD format. This format allows for formal charges, implicit hydrogen, bond orders and total molecular charges to be defined. It has three possible entry types:
- "ATOM  X  FC  H". This type for specifing things to do with atoms. The header 'ATOM' must always be included. X is the atom index which the settings should be applied to, FC is the formal charge and H is the number of implicit hydrogens.
- "BOND  X  Y  O". This type is for specifing bond order between two atoms. The header 'BOND' must always be included. X and Y are the atom indices and O is the order. O can be 1 - 4 for single to quadruple bonds, 5 for aromatic, 6 for 1.5 and 7 for 2.5.
- "MOLECULE  Q". This type is for specifing the total molecular charge of a molecule. The header 'MOLECULE' must always be included. Q is the total molecular charge. This must match the sum of formal charges otherwise an error will be thrown.

The second file format is the fragment definition file. This format contains three blocks, each opened with a header and terminated with 'END'.
- The 'MOLECULE' block must always come first in the file, and must only occur once. It contains three lines indicating which files to load for the molecule definition. First, a coordinate file, second a parameter file and last an IXD file. These file paths are relative to the fragment file path.
- The 'FRAGMENT' block contains the list of atom indices to use for the core region in fragment generation. This is a free-form block taking any number of lines with any number of indices per line.
- The 'OVERLAP' block contains the list of atom inidices to use for the overlap region in fragment generation. Like the 'FRAGMENT' block, this is free form. Each occurance of a 'FRAGMENT' block requires the next block to be an 'OVERLAP' block, even if the 'OVERLAP' block is empty.

