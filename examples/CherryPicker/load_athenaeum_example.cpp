//
// Created by sdun067 on 12/11/19.
//

#include <indigox/indigox.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/forcefield.hpp>


int main() {

  using namespace indigox;

  auto forceField = GenerateGROMOS54A7();
//  auto athSettings = Athenaeum::Settings();

  auto auto_ath = Athenaeum(forceField, 1);
  auto_ath.SetInt(Athenaeum::Settings::MoleculeSizeLimit, 60);

  //construct molecule. Simple is fine.
  
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

  auto_ath.AddAllFragments(mol);

  SaveAthenaeum(auto_ath, "attemptedLib.ath");

  std::cout << "Done!";
}
