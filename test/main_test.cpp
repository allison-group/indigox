//#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
//#include <indigox/utils/doctest_proxy.hpp>
//
//#include <indigox/utils/common.hpp>
//#include <indigox/utils/counter.hpp>
//
//#include <indigox/algorithm/graph/connectivity.hpp>
//#include <indigox/algorithm/graph/paths.hpp>
//#include <indigox/classes/athenaeum.hpp>
//#include <indigox/classes/molecule.hpp>
//#include <indigox/graph/molecular.hpp>
//#include <indigox/graph/condensed.hpp>
//
////test_case("IXCountableObject counts correctly") {
////  auto int_a = indigox::utils::IXCountableObject<int>();
////  auto int_b = indigox::utils::IXCountableObject<int>();
////  auto float_a = indigox::utils::IXCountableObject<float>();
////  auto float_b = indigox::utils::IXCountableObject<float>();
////
////  check(int_a.GetUniqueID() + 1 == int_b.GetUniqueID());
////  check(float_a.GetUniqueID() + 1 == float_b.GetUniqueID());
////  check(int_a.GetUniqueID() == float_a.GetUniqueID());
////}
//
//
//test_case("ShortTest") {
//  using namespace indigox;
//  Molecule mol = CreateMolecule();
//  Atom C1 = mol->NewAtom("C1", GetPeriodicTable()->GetElement("C"));
//  Atom C2 = mol->NewAtom("C2", GetPeriodicTable()->GetElement("C"));
//  Atom C3 = mol->NewAtom("C3", GetPeriodicTable()->GetElement("C"));
//  Atom C4 = mol->NewAtom("C4", GetPeriodicTable()->GetElement("C"));
//  Atom C5 = mol->NewAtom("C5", GetPeriodicTable()->GetElement("C"));
//  Atom C6 = mol->NewAtom("C6", GetPeriodicTable()->GetElement("C"));
//  Atom H1 = mol->NewAtom("H1", GetPeriodicTable()->GetElement("H"));
//  Atom Cl = mol->NewAtom("Cl", GetPeriodicTable()->GetElement("Cl"));
//  Atom F1 = mol->NewAtom("F1", GetPeriodicTable()->GetElement("F"));
//  Atom Br = mol->NewAtom("Br", GetPeriodicTable()->GetElement("Br"));
//  Atom I1 = mol->NewAtom("I1", GetPeriodicTable()->GetElement("I"));
//  Atom H2 = mol->NewAtom("H2", GetPeriodicTable()->GetElement("H"));
//  mol->NewBond(C1, C2)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C2, C5)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C2, C3)->SetOrder(BondOrder::DOUBLE);
//  mol->NewBond(C3, C4)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C4, C5)->SetOrder(BondOrder::DOUBLE);
//  mol->NewBond(C5, C6)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C6, C1)->SetOrder(BondOrder::DOUBLE);
//  mol->NewBond(C1, H1)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C2, Cl)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C3, F1)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C4, Br)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C5, I1)->SetOrder(BondOrder::SINGLE);
//  mol->NewBond(C6, H2)->SetOrder(BondOrder::SINGLE);
//
//  Molecule mol2 = CreateMolecule();
//  Atom C11 = mol2->NewAtom("C1", GetPeriodicTable()->GetElement("C"));
//  Atom C21 = mol2->NewAtom("C2", GetPeriodicTable()->GetElement("C"));
//  Atom C31 = mol2->NewAtom("C3", GetPeriodicTable()->GetElement("C"));
//  Atom C41 = mol2->NewAtom("C4", GetPeriodicTable()->GetElement("C"));
//  Atom C51 = mol2->NewAtom("C5", GetPeriodicTable()->GetElement("C"));
//  Atom C61 = mol2->NewAtom("C6", GetPeriodicTable()->GetElement("C"));
//  Atom H11 = mol2->NewAtom("H1", GetPeriodicTable()->GetElement("H"));
//  Atom Cl1 = mol2->NewAtom("Cl", GetPeriodicTable()->GetElement("Cl"));
//  Atom F11 = mol2->NewAtom("F1", GetPeriodicTable()->GetElement("F"));
//  Atom Br1 = mol2->NewAtom("Br", GetPeriodicTable()->GetElement("Br"));
//  Atom I11 = mol2->NewAtom("I1", GetPeriodicTable()->GetElement("I"));
//  Atom H21 = mol2->NewAtom("H2", GetPeriodicTable()->GetElement("H"));
//  mol2->NewBond(C11, C21)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C21, C31)->SetOrder(BondOrder::DOUBLE);
//  mol2->NewBond(C31, C41)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C41, C51)->SetOrder(BondOrder::DOUBLE);
//  mol2->NewBond(C51, C61)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C61, C11)->SetOrder(BondOrder::DOUBLE);
//  mol2->NewBond(C11, H11)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C21, Cl1)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C31, F11)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C41, Br1)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C51, I11)->SetOrder(BondOrder::SINGLE);
//  mol2->NewBond(C61, H21)->SetOrder(BondOrder::SINGLE);
//
//  Athenaeum ath = std::make_shared<IXAthenaeum>(2, 1);
//  std::cout << ath->AddAllFragments(mol) << "\n";
//  std::cout << ath->AddAllFragments(mol) << "\n";
//  std::cout << ath->NumFragments() << "\n";
//  std::cout << ath->NumFragments(mol) << "\n";
//  std::cout << ath->AddAllFragments(mol2) << "\n";
//  std::cout << ath->AddAllFragments(mol) << "\n";
//  std::cout << ath->NumFragments() << "\n";
//  std::cout << ath->NumFragments(mol) << "\n";
//  std::cout << ath->NumFragments(mol2) << "\n";
//
//}
#include <iostream>
#include <memory>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/molecular.hpp>

struct A {
  int dat;
  A(int i) : dat(i) { }
};

struct B {
  A& a;
  B(A& _a) : a(_a) { }
};

A& makeA() { return *std::make_shared<A>(5); }

int main() {
  using namespace indigox;
  
  sMolecule mol = CreateMolecule();
  Atom& a1 = mol->NewAtom("A1");
  Atom& a2 = mol->NewAtom("A2");
  Atom& a3 = mol->NewAtom("A3");
  graph::MolecularGraph& G = mol->GetGraph();
//  std::cout << "Connected: " << G.NumConnectedComponents() << ", Cycles: " << G.NumCycles() << "\n";
  Bond& b1 = mol->NewBond(a1, a2);
//  std::cout << "Connected: " << G.NumConnectedComponents() << ", Cycles: " << G.NumCycles() << "\n";
  Bond& b2 = mol->NewBond(a1, a3);
  std::cout << "Connected: " << G.NumConnectedComponents() << ", Cycles: " << G.NumCycles() << "\n";
  Bond& b3 = mol->NewBond(a3, a2);
  std::cout << "Connected: " << G.NumConnectedComponents() << ", Cycles: " << G.NumCycles() << "\n";
  
  return 0;
}
