#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>

#include <indigox/utils/common.hpp>

#include "class_test_wrappers.hpp"

namespace indigox::test {
  struct MoleculeFixture {
    indigox::test::IXMolecule methane, mol;
    PeriodicTable PT;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    MoleculeFixture() : methane(), mol(), PT(GetPeriodicTable()) {
      // Build methane
      atoms.emplace_back(methane.NewAtom());
      atoms.back()->SetElement("C");
      atoms.back()->SetTag(0);
      for (size_ i = 0; i < 4; ++i) {
        atoms.emplace_back(methane.NewAtom());
        if (i < 2) atoms.back()->SetElement("F");
        else if (i == 2) atoms.back()->SetElement("H");
        else atoms.back()->SetElement("Cl");
        bonds.emplace_back(methane.NewBond(atoms.front(), atoms.back()));
        // Tags will be: 1, 2, 5, 10
        atoms.back()->SetTag(i * i + 1);
        bonds.back()->SetTag(i * i + 1);
      }
      bonds.front()->SwapSourceTarget();  // Help catch all regions of HasBond
    }
  };
}

using namespace indigox;

BOOST_FIXTURE_TEST_SUITE(ixmolecule, indigox::test::MoleculeFixture);

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK_NO_THROW(indigox::test::IXMolecule());
  
  // All emergent properties should be calculated
  BOOST_CHECK(methane.GetEmergentState().all());
  
  BOOST_CHECK(methane.NumAtoms() == 5);
  BOOST_CHECK(methane.NumBonds() == 4);
//  BOOST_CHECK(methane.NumAngles() == 0);
//  BOOST_CHECK(methane.NumDihedrals() == 0);
  
  methane.ReserveAtoms(200);
  BOOST_CHECK(methane.AtomCapacity() >= 200);
  methane.ReserveBonds(2000);
  BOOST_CHECK(methane.BondCapacity() >= 2000);
}

BOOST_AUTO_TEST_CASE(name_get_set) {
  BOOST_CHECK(mol.GetName() == "");
  mol.SetName("Testerowni");
  BOOST_CHECK(mol.GetName() == "Testerowni");
}

BOOST_AUTO_TEST_CASE(molecular_charge_get_set) {
  methane.ResetEmergentState();
  BOOST_CHECK(methane.GetMolecularCharge() == 0);
  test::IXMolecule::EmergeSet state = methane.GetEmergentState();
  methane.SetMolecularCharge(12); // Should change emergent state
  BOOST_CHECK(state != methane.GetEmergentState());
  BOOST_CHECK(methane.GetMolecularCharge() == 12);
  methane.ResetEmergentState();
  // Setting to current value should not change emergent state
  methane.SetMolecularCharge(12);
  BOOST_CHECK(state == methane.GetEmergentState());
}

BOOST_AUTO_TEST_CASE(atom_has) {
  std::vector<bool> expected, obtained;
  for (Atom atm : atoms) {
    expected.push_back(true); expected.push_back(false);
    obtained.push_back(methane.HasAtom(atm));
    obtained.push_back(mol.HasAtom(atm));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                obtained.begin(), obtained.end());
  BOOST_CHECK(!mol.HasAtom(Atom()));
}

BOOST_AUTO_TEST_CASE(bond_has) {
  // Check all expected bonds exist and no bonus ones, via atom pairs
  std::vector<bool> expected = {true, true, true, true, false, false, false,
                                false, false, false};
  std::vector<bool> obtained;
  for (size_ i = 0; i < atoms.size(); ++i) {
    for (size_ j = i + 1; j < atoms.size(); ++ j)
      obtained.push_back(methane.HasBond(atoms[j], atoms[i]));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                obtained.begin(), obtained.end());
  
  // Check bonds via bonds
  expected.clear(); obtained.clear();
  for (Bond bnd : bonds) {
    expected.push_back(true); expected.push_back(false);
    obtained.push_back(methane.HasBond(bnd));
    obtained.push_back(mol.HasBond(bnd));
  }
  
  BOOST_CHECK(!mol.HasBond(Bond()));
  BOOST_CHECK(!mol.HasBond(Atom(), Atom()));
  BOOST_CHECK(!mol.HasBond(atoms[0], atoms[1]));
}


BOOST_AUTO_TEST_CASE(atom_get) {
  std::vector<Atom> got_atoms;
  // Get by position
  for (size_ p = 0; p < methane.NumAtoms(); ++p)
    got_atoms.emplace_back(methane.GetAtom(p));
  BOOST_CHECK_EQUAL_COLLECTIONS(atoms.begin(), atoms.end(),
                                got_atoms.begin(), got_atoms.end());
  got_atoms.clear();
  // over position should return null
  BOOST_CHECK(methane.GetAtom(methane.NumAtoms()) == Atom());
  
  
  // Get by id
  for (size_ p = 0; p < methane.NumAtoms(); ++p)
    got_atoms.emplace_back(methane.GetAtomID(atoms[p]->GetUniqueID()));
  BOOST_CHECK_EQUAL_COLLECTIONS(atoms.begin(), atoms.end(),
                                got_atoms.begin(), got_atoms.end());
  got_atoms.clear();
  // Bad id should return null
  BOOST_CHECK(methane.GetAtomID(std::numeric_limits<uid_>::max()) == Atom());
  
  // Get by tag
  for (size_ p = 0; p < methane.NumAtoms(); ++p)
    got_atoms.emplace_back(methane.GetAtomTag(atoms[p]->GetTag()));
  BOOST_CHECK_EQUAL_COLLECTIONS(atoms.begin(), atoms.end(),
                                got_atoms.begin(), got_atoms.end());
  got_atoms.clear();
  // Get by tag should get the first atom with the given tag
  atoms.front()->SetTag(atoms.back()->GetTag());
  BOOST_CHECK(methane.GetAtomTag(atoms.back()->GetTag()) == atoms.front());
  // Bad tag should return null
  BOOST_CHECK(methane.GetAtomTag(900) == Atom());
  
  // Get by iterator
  auto atm_itr = methane.GetAtoms();
  BOOST_CHECK_EQUAL_COLLECTIONS(atm_itr.first, atm_itr.second,
                                atoms.begin(), atoms.end());
}

BOOST_AUTO_TEST_CASE(bond_get) {
  std::vector<Bond> got_bonds;
  // Get by position
  for (size_ p = 0; p < methane.NumBonds(); ++p)
    got_bonds.emplace_back(methane.GetBond(p));
  BOOST_CHECK_EQUAL_COLLECTIONS(got_bonds.begin(), got_bonds.end(),
                                bonds.begin(), bonds.end());
  got_bonds.clear();
  // over position should return null
  BOOST_CHECK(methane.GetBond(methane.NumBonds()) == Bond());
  
  // Get by id
  for (size_ p = 0; p < methane.NumBonds(); ++p)
    got_bonds.emplace_back(methane.GetBondID(bonds[p]->GetUniqueID()));
  BOOST_CHECK_EQUAL_COLLECTIONS(got_bonds.begin(), got_bonds.end(),
                                bonds.begin(), bonds.end());
  got_bonds.clear();
  // Bad id should return null
  BOOST_CHECK(methane.GetBondID(std::numeric_limits<uid_>::max()) == Bond());
  
  // Get by tag
  for (size_ p = 0; p < methane.NumBonds(); ++p)
    got_bonds.emplace_back(methane.GetBondTag(bonds[p]->GetTag()));
  BOOST_CHECK_EQUAL_COLLECTIONS(got_bonds.begin(), got_bonds.end(),
                                bonds.begin(), bonds.end());
  got_bonds.clear();
  // Get by tag should get the first bond with the given tag
  bonds.front()->SetTag(bonds.back()->GetTag());
  BOOST_CHECK(methane.GetBondTag(bonds.back()->GetTag()) == bonds.front());
  // Bad tag should return null
  BOOST_CHECK(methane.GetBondTag(900) == Bond());
  
  // Get by iterator
  auto bnd_itr = methane.GetBonds();
  BOOST_CHECK_EQUAL_COLLECTIONS(bnd_itr.first, bnd_itr.second,
                                bonds.begin(), bonds.end());
  
  // Get by atom pairs
  std::vector<Bond> expected = {bonds[0], bonds[1], bonds[2], bonds[3], Bond(),
                                Bond(), Bond(), Bond(), Bond(), Bond()};
  for (size_ i = 0; i < atoms.size(); ++i) {
    for (size_ j = i + 1; j < atoms.size(); ++ j)
      got_bonds.push_back(methane.GetBond(atoms[j], atoms[i]));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                got_bonds.begin(), got_bonds.end());
  
}

BOOST_AUTO_TEST_CASE(atom_new_remove) {
  // Normal new atom
  Atom atm = mol.NewAtom();
  BOOST_CHECK(atm);
  BOOST_CHECK(mol.NumAtoms() == 1);
  BOOST_CHECK(utils::IsNull(atm->GetElement()));
  
  // New atom with element
  atm = mol.NewAtom(PT->GetElement(6));
  BOOST_CHECK(atm);
  BOOST_CHECK(mol.NumAtoms() == 2);
  BOOST_CHECK(PT->GetElement(6) == atm->GetElement());
  
  // New atom with name
  atm = mol.NewAtom("Test");
  BOOST_CHECK(atm);
  BOOST_CHECK(mol.NumAtoms() == 3);
  BOOST_CHECK(atm->GetName() == "Test");
  
  // New atom with name and element
  atm = mol.NewAtom("Test2", PT->GetElement(12));
  BOOST_CHECK(atm);
  BOOST_CHECK(mol.NumAtoms() == 4);
  BOOST_CHECK(atm->GetName() == "Test2");
  BOOST_CHECK(atm->GetElement() == PT->GetElement(12));
  
  // Remove atom
  BOOST_CHECK(mol.RemoveAtom(atm));
  BOOST_CHECK(mol.NumAtoms() == 3);
  
  // Remove already removed should fail
  BOOST_CHECK(!mol.RemoveAtom(atm));
  BOOST_CHECK(mol.NumAtoms() == 3);
  
  // Remove not part of should fail
  BOOST_CHECK(!mol.RemoveAtom(atoms.front()));
  BOOST_CHECK(mol.NumAtoms() == 3);
  
  // Remove null should fail
  BOOST_CHECK(!mol.RemoveAtom(Atom()));
  BOOST_CHECK(mol.NumAtoms() == 3);
}

BOOST_AUTO_TEST_CASE(bond_new_remove) {
  Atom atm1 = mol.NewAtom();
  Atom atm2 = mol.NewAtom();
  Bond bnd1 = mol.NewBond(atm1, atm2);
  BOOST_CHECK(bnd1);
  BOOST_CHECK(mol.NumBonds() == 1);
  BOOST_CHECK(bnd1->GetSourceAtom() == atm1);
  BOOST_CHECK(bnd1->GetTargetAtom() == atm2);
  
  // Adding bond between atoms that already have a bonds should return null_ptr
  Bond bnd2 = mol.NewBond(atm2, atm1);
  BOOST_CHECK(!bnd2);
  BOOST_CHECK(mol.NumBonds() == 1);
  
  // Adding bond with not owned atom should return null_ptr
  Bond bnd3 = mol.NewBond(atoms.front(), atoms.back());
  BOOST_CHECK(!bnd3);
  BOOST_CHECK(mol.NumBonds() == 1);
  
  mol.NewBond(mol.NewAtom(), mol.NewAtom());
  
  // Remove bond
  BOOST_CHECK(mol.RemoveBond(bnd1));
  BOOST_CHECK(mol.NumBonds() == 1);
  
  // Remove already removed should fail
  BOOST_CHECK(!mol.RemoveBond(bnd1));
  BOOST_CHECK(mol.NumBonds() == 1);
  
  // Remove not part of should fail
  BOOST_CHECK(!mol.RemoveBond(bonds.front()));
  BOOST_CHECK(mol.NumBonds() == 1);
  
  // Remove null should fail
  BOOST_CHECK(!mol.RemoveBond(Bond()));
  BOOST_CHECK(mol.NumBonds() == 1);
  
  // Remove bond from between atoms
  BOOST_CHECK(methane.RemoveBond(atoms.front(), atoms.back()));
  BOOST_CHECK(methane.NumBonds() == 3);
  
  // Remove already removed bond should fail
  BOOST_CHECK(!methane.RemoveBond(atoms.front(), atoms.back()));
  BOOST_CHECK(methane.NumBonds() == 3);
  
  // Remove never existed bond should fail
  BOOST_CHECK(!methane.RemoveBond(atoms[2], atoms[1]));
  BOOST_CHECK(methane.NumBonds() == 3);
  
  // Remove not owned atoms should fail
  BOOST_CHECK(!mol.RemoveBond(atoms[1], atoms.front()));
  BOOST_CHECK(mol.NumBonds() == 1);
  
  // Remove null atoms should fail
  BOOST_CHECK(!mol.RemoveBond(Atom(), Atom()));
  BOOST_CHECK(mol.NumBonds() == 1);
  
}

BOOST_AUTO_TEST_CASE(formula_get) {
  using Emergent = indigox::test::IXMolecule::Emergent;
  size_ pos = static_cast<size_>(Emergent::MOLECULAR_FORMULA);
  
  // Check correct calculation and reset Emergent::MOLECULAR_FORMULA on call
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  BOOST_CHECK(methane.GetFormula() == "CHClF2");
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  atoms.front()->SetElement("U");
  BOOST_CHECK(methane.GetFormulaCache() == "CHClF2");
  BOOST_CHECK(methane.GetFormula() == "HClF2U");
  BOOST_CHECK(methane.GetFormulaCache() == "HClF2U");
  for (Atom atm : atoms) atm->SetElement("C");
  BOOST_CHECK(methane.GetFormula() == "C5");
  for (Atom atm : atoms) atm->SetElement("H");
  BOOST_CHECK(methane.GetFormula() == "H5");
  
  // Check only calculates when Emergent::MOLECULAR_FORMULA is set
  atoms.front()->SetElement("Si");
  methane.ResetEmergentState(pos); // should stop recalc
  BOOST_CHECK(methane.GetFormula() == "H5"); // Should return cached
  methane.SetEmergentState(pos); // should enable recalc
  BOOST_CHECK(methane.GetFormula() == "H4Si"); // Should recalc and return
}

BOOST_AUTO_TEST_CASE(graph_get_check) {
  graph::MolecularGraph G = methane.GetGraph();
  BOOST_CHECK(G != graph::MolecularGraph());
  BOOST_CHECK(G->NumEdges() == methane.NumBonds());
  BOOST_CHECK(G->NumVertices() == methane.NumAtoms());
  
  // Vertices should reference all atoms, edges reference all bonds
  auto atm_itr = methane.GetAtoms(); auto bnd_itr = methane.GetBonds();
  auto vert_it = G->GetVertices(); auto edge_it = G->GetEdges();
  std::set<Atom> expected_at(atm_itr.first, atm_itr.second);
  std::set<Bond> expected_bn(bnd_itr.first, bnd_itr.second);
  std::set<Atom> obtained_at;
  std::set<Bond> obtained_bn;
  for (; vert_it.first != vert_it.second; ++vert_it.first)
    obtained_at.emplace((*vert_it.first)->GetAtom());
  for (; edge_it.first != edge_it.second; ++edge_it.first)
    obtained_bn.emplace((*edge_it.first)->GetBond());
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_at.begin(), expected_at.end(),
                                obtained_at.begin(), obtained_at.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_bn.begin(), expected_bn.end(),
                                obtained_bn.begin(), obtained_bn.end());
  
  // Adding an atom should add a vertex
  Atom test1 = methane.NewAtom();
  BOOST_CHECK(G->NumVertices() == methane.NumAtoms());
  
  // Adding a bond should add an edge
  Bond test2 = methane.NewBond(test1, atoms.back());
  BOOST_CHECK(G->NumEdges() == methane.NumBonds());
  
  // Removing a bond should remove an edge
  methane.RemoveBond(test2);
  BOOST_CHECK(G->NumEdges() == methane.NumBonds());
  
  // Removing a bond between atoms should remove an edge
  methane.RemoveBond(atoms.back(), atoms.front());
  BOOST_CHECK(G->NumEdges() == methane.NumBonds());
  
  // Removing an atom should remove a vertex
  methane.RemoveAtom(test1);
  BOOST_CHECK(G->NumVertices() == methane.NumAtoms());
  
  // Removing an atom with bonds should remove vertex and edges
  methane.RemoveAtom(atoms.front());
  BOOST_CHECK(G->NumVertices() == methane.NumAtoms());
  BOOST_CHECK(G->NumEdges() == methane.NumBonds());
}

BOOST_AUTO_TEST_SUITE(emergent);

BOOST_AUTO_TEST_CASE(molecular_formula) {
  using Emergent = indigox::test::IXMolecule::Emergent;
  size_ pos = static_cast<size_>(Emergent::MOLECULAR_FORMULA);
  methane.ResetEmergentState();
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  
  // Adding an atom should set
  Atom test1 = methane.NewAtom();
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState();
  
  // Removing an atom should set
  methane.RemoveAtom(test1);
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState();
  
  // Changing an atom's element should set
  atoms.front()->SetElement("Se");
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState();
}

BOOST_AUTO_TEST_CASE(topological_bofc) {
  using Emergent = indigox::test::IXMolecule::Emergent;
  size_ pos = static_cast<size_>(Emergent::TOPOLOGICAL_BOFC);
  methane.ResetEmergentState(pos);
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  
  // Adding an atom should set
  Atom test1 = methane.NewAtom();
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Adding a bond should set
  Bond test2 = methane.NewBond(test1, atoms.front());
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing a bond should set
  methane.RemoveBond(test2);
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing a bond between atoms should also set
  methane.RemoveBond(atoms.front(), atoms.back());
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing an atom should set
  methane.RemoveAtom(test1);
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Changing an element should set
  atoms.front()->SetElement("Xe");
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Changing the molecular charge should set
  methane.SetMolecularCharge(25);
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
}

BOOST_AUTO_TEST_CASE(angle_perception) {
  using Emergent = indigox::test::IXMolecule::Emergent;
  size_ pos = static_cast<size_>(Emergent::ANGLE_PERCEPTION);
  methane.ResetEmergentState(pos);
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  
  // Adding an atom should not set
  Atom test1 = methane.NewAtom();
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Adding a bond should set
  Bond test2 = methane.NewBond(test1, atoms.front());
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing a bond should set
  methane.RemoveBond(test2);
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing a bond between atoms should also set
    methane.RemoveBond(atoms.front(), atoms.back());
    BOOST_CHECK(methane.GetEmergentState()[pos]);
    methane.ResetEmergentState(pos);
  
  // Removing an atom should not set
  methane.RemoveAtom(test1);
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing an atom with bonds should set
  methane.RemoveAtom(atoms.front());
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
}

BOOST_AUTO_TEST_CASE(dihedral_perception) {
  using Emergent = indigox::test::IXMolecule::Emergent;
  size_ pos = static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION);
  methane.ResetEmergentState(pos);
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  
  // Adding an atom should not set
  Atom test1 = methane.NewAtom();
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Adding a bond should set
  Bond test2 = methane.NewBond(test1, atoms.front());
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing a bond should set
  methane.RemoveBond(test2);
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing a bond between atoms should also set
    methane.RemoveBond(atoms.front(), atoms.back());
    BOOST_CHECK(methane.GetEmergentState()[pos]);
    methane.ResetEmergentState(pos);
  
  // Removing an atom should not set
  methane.RemoveAtom(test1);
  BOOST_CHECK(!methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
  
  // Removing an atom with bonds should set
  methane.RemoveAtom(atoms.front());
  BOOST_CHECK(methane.GetEmergentState()[pos]);
  methane.ResetEmergentState(pos);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_CASE(property_modified) {
  using Emergent = indigox::test::IXMolecule::Emergent;
  size_ form = static_cast<size_>(Emergent::MOLECULAR_FORMULA);
  size_ topo = static_cast<size_>(Emergent::TOPOLOGICAL_BOFC);
  size_ angp = static_cast<size_>(Emergent::ANGLE_PERCEPTION);
  size_ dhdp = static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION);
  mol.ResetEmergentState();
  
  // ATOM_ELEMENTS should set MOLECULAR_FORMULA, TOPOLOGICL_BOFC
  mol.SetPropertyModified(MolProperty::ATOM_ELEMENTS);
  BOOST_CHECK(mol.GetEmergentState()[form]);
  BOOST_CHECK(mol.GetEmergentState()[topo]);
  BOOST_CHECK(!mol.GetEmergentState()[angp]);
  BOOST_CHECK(!mol.GetEmergentState()[dhdp]);
  mol.ResetEmergentState();
  
  // CONNECTIVITY should set TOPOLOGICL_BOFC, angle/dihed perception
  mol.SetPropertyModified(MolProperty::CONNECTIVITY);
  BOOST_CHECK(!mol.GetEmergentState()[form]);
  BOOST_CHECK(mol.GetEmergentState()[topo]);
  BOOST_CHECK(mol.GetEmergentState()[angp]);
  BOOST_CHECK(mol.GetEmergentState()[dhdp]);
  mol.ResetEmergentState();
  
  // ELECTRON_COUNT should set TOPOLOGICL_BOFC
  mol.SetPropertyModified(MolProperty::ELECTRON_COUNT);
  BOOST_CHECK(!mol.GetEmergentState()[form]);
  BOOST_CHECK(mol.GetEmergentState()[topo]);
  BOOST_CHECK(!mol.GetEmergentState()[angp]);
  BOOST_CHECK(!mol.GetEmergentState()[dhdp]);
  mol.ResetEmergentState();
  
  // NUM_PROPERTIES should set NUM_EMERGENT (and never be used...)
  mol.SetPropertyModified(MolProperty::NUM_PROPERTIES);
  BOOST_CHECK(!mol.GetEmergentState()[form]);
  BOOST_CHECK(!mol.GetEmergentState()[topo]);
  BOOST_CHECK(!mol.GetEmergentState()[angp]);
  BOOST_CHECK(!mol.GetEmergentState()[dhdp]);
}


BOOST_AUTO_TEST_SUITE_END();
