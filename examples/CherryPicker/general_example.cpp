//
// C++ example that imports a molecule from an outputted binary (which you can produce from python with indigox.SaveMolecule)
// This example has CalculateElectrons on, so it will calculate formal charges and bond orders.
// Relies on the Python example having been run first so there is an Athenaeum to match the molecule to.
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

  auto forceField = GenerateGROMOS54A7();

  auto man_ath = LoadAthenaeum(example_folder_path + "ManualAthenaeum.ath");
//  auto auto_ath = LoadAthenaeum(example_folder_path + "AutomaticAthenaeum.ath");
  Molecule mol = LoadMolecule(example_folder_path + "TestMolecules/axinellinA.bin");

  algorithm::CherryPicker cherryPicker(forceField);
  cherryPicker.AddAthenaeum(man_ath);

  cherryPicker.SetInt(settings::MinimumFragmentSize, 2);
  cherryPicker.SetInt(settings::MaximumFragmentSize, 20);
  cherryPicker.SetBool(settings::CalculateElectrons);

  const indigox::ParamMolecule &molecule = cherryPicker.ParameteriseMolecule(mol);

  //We can't save the parameterisation directly in ITP format from C++, but we save the binary molecule output
  //which can be imported by the python module and outputted in ITP, IXD, PDB and RTP formats
  SaveMolecule(mol, example_folder_path + "TestMolecules/axinellinA.out.param");

  std::cout << "Done!\n";
}
