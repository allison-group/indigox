//
// Pure C++ example, useful for if you need to debug small issues in C++.
// Relies on the Python example having been run first so there is an Athenaeum to match the 4-Nitrobenzoate molecule to.
// Also assumes you are running from a build folder in the project root. If not, change example_folder_path below
//

#include <indigox/indigox.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/parameterised.hpp>
#include <indigox/algorithm/cherrypicker.hpp>
#include <experimental/filesystem>


int main() {

  using namespace indigox;
  namespace fs = std::experimental::filesystem;
  using settings = indigox::algorithm::CherryPicker::Settings;

  std::string example_folder_path = "../examples/CherryPicker/";

  // Always need a forcefield. The molecule, Athenaeum and CherryPicker forcefields must all match
  auto forceField = GenerateGROMOS54A7();

  //construct a simple molecule

  // Prepare elements
  const PeriodicTable& PT = GetPeriodicTable();
  Element H = PT.GetElement("H");
  Element C = PT.GetElement("C");
  Element O = PT.GetElement("O");
  Element N = PT.GetElement("N");

  Molecule mol = Molecule("4-Nitrobenzoate");
  mol.SetForcefield(forceField);

  Atom c1  = mol.NewAtom(C);  c1.SetName("C1") ;
  Atom c2  = mol.NewAtom(C);  c2.SetName("C2") ;
  Atom c3  = mol.NewAtom(C);  c3.SetName("C3") ;
  Atom c4  = mol.NewAtom(C);  c4.SetName("C4") ;
  Atom c5  = mol.NewAtom(C);  c5.SetName("C5") ;
  Atom c6  = mol.NewAtom(C);  c6.SetName("C6") ;
  Atom h7  = mol.NewAtom(H);  h7.SetName("H7") ;
  Atom h8  = mol.NewAtom(H);  h8.SetName("H8") ;
  Atom n9  = mol.NewAtom(N);  n9.SetName("N9") ;
  Atom h10 = mol.NewAtom(H); h10.SetName("H10");
  Atom h11 = mol.NewAtom(H); h11.SetName("H11");
  Atom c12 = mol.NewAtom(C); c12.SetName("C12");
  Atom o13 = mol.NewAtom(O); o13.SetName("O13");
  Atom o14 = mol.NewAtom(O); o14.SetName("O14");
  Atom o15 = mol.NewAtom(O); o15.SetName("O15");
  Atom o16 = mol.NewAtom(O); o16.SetName("O16");

  mol.NewBond(c1,c2);
  mol.NewBond(c1,c6);
  mol.NewBond(c1,h7);
  mol.NewBond(c2,c3);
  mol.NewBond(c2,h8);
  mol.NewBond(c3,c4);
  mol.NewBond(c3,n9);
  mol.NewBond(c4,c5);
  mol.NewBond(c4,h10);
  mol.NewBond(c5,c6);
  mol.NewBond(c5,h11);
  mol.NewBond(c6,c12);
  mol.NewBond(n9,o13);
  mol.NewBond(n9,o14);
  mol.NewBond(c12,o15);
  mol.NewBond(c12,o16);

  // Load the Athenaeum
  auto auto_ath = LoadAthenaeum(example_folder_path + "AutomaticAthenaeum.ath");

  // Set up CherryPicker
  algorithm::CherryPicker cherryPicker(forceField);
  cherryPicker.AddAthenaeum(auto_ath);
  
  cherryPicker.SetInt(settings::MinimumFragmentSize, 2);
  cherryPicker.SetInt(settings::MaximumFragmentSize, 20);

  const indigox::ParamMolecule &molecule = cherryPicker.ParameteriseMolecule(mol);

  SaveMolecule(mol, "../examples/CherryPicker/TestMolecules/nitrobenzoate.out");

  std::cout << "Done!\n";
}
