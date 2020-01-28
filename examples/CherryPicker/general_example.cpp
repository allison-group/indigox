//
// Created by sdun067 on 12/11/19.
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

  auto forceField = GenerateGROMOS54A7();
  auto athSettings = Athenaeum::Settings();

  auto man_ath = LoadAthenaeum("/home/sdun067/AllisonGroup/indigox/examples/CherryPicker/ManualAthenaeum.ath");
//  auto auto_ath = LoadAthenaeum("/home/sdun067/AllisonGroup/indigox/examples/CherryPicker/AutomaticAthenaeum.ath");
  Molecule mol = LoadMolecule("/home/sdun067/AllisonGroup/indigox/examples/CherryPicker/TestMolecules/polymyxin.out");
//  Molecule mol = LoadMolecule("/home/sdun067/AllisonGroup/indigox/examples/CherryPicker/TestMolecules/nitrobenzoate.out");

  algorithm::CherryPicker cherryPicker(forceField);
  cherryPicker.AddAthenaeum(man_ath);
//  cherryPicker.AddAthenaeum(auto_ath);

  cherryPicker.SetInt(settings::MinimumFragmentSize, 2);
  cherryPicker.SetInt(settings::MaximumFragmentSize, 20);
  cherryPicker.SetInt(settings::ElectronMethod, 2);

  const indigox::ParamMolecule &molecule = cherryPicker.ParameteriseMolecule(mol);

  SaveMolecule(mol, "/home/sdun067/AllisonGroup/indigox/examples/CherryPicker/TestMolecules/general.out.param");

  std::cout << "Done!\n";
}
