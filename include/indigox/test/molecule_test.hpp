#ifndef INDIGOX_TEST_MOLECULE_TEST_HPP
#define INDIGOX_TEST_MOLECULE_TEST_HPP

#include <random>

namespace indigox::test {
  struct TestMolecule {
    // Typedefs
    using IAtoms = IXMolecule::MolAtoms;
    using IBonds = IXMolecule::MolBonds;
    using IAngles = IXMolecule::MolAngles;
    using IDihedrals = IXMolecule::MolDihedrals;
    using Atoms = std::pair<IXMolecule::MolAtomIter, IXMolecule::MolAtomIter>;
    using Bonds = std::pair<IXMolecule::MolBondIter, IXMolecule::MolBondIter>;
    using Angles = std::pair<IXMolecule::MolAngleIter, IXMolecule::MolAngleIter>;
    using Dihedrals = std::pair<IXMolecule::MolDihedralIter, IXMolecule::MolDihedralIter>;
    using EmergeSet = IXMolecule::EmergeSet;
    using EmergeProp = IXMolecule::Emergent;
    
    indigox::Molecule imp;
  
    // Private wrapping functions
    TestMolecule() : imp(new IXMolecule()) { }
    TestMolecule(bool) : imp(CreateMolecule()) { }
    void Init() { imp->Init(); }
    Angle NewAngle(const Atom& a, const Atom& b, const Atom& c) {
      return imp->NewAngle(a,b,c);
    }
    Dihedral NewDihedral(const Atom& a, const Atom& b,
                         const Atom& c, const Atom& d) {
      return imp->NewDihedral(a,b,c,d);
    }
    Bond FindBond(const Atom& a, const Atom& b) const {
      return imp->_FindBond(a,b);
    }
    Angle FindAngle(const Atom& a, const Atom& b, const Atom& c) const {
      return imp->_FindAngle(a,b,c);
    }
    Dihedral FindDihedral(const Atom& a, const Atom& b,
                          const Atom& c, const Atom& d) const {
      return imp->_FindDihedral(a,b,c,d);
    }
    
    // Public wrapping functions
    Angle GetAngle(size_ p) const { return imp->GetAngle(p); }
    Dihedral GetDihedral(size_ p) const { return imp->GetDihedral(p); }
    Angle GetAngle(const Atom& a, const Atom& b, const Atom& c) {
      return imp->GetAngle(a,b,c);
    }
    Dihedral GetDihedral(const Atom& a, const Atom& b,
                         const Atom& c, const Atom& d) {
      return imp->GetDihedral(a,b,c,d);
    }
    Angle GetAngleTag(uid_ t) const { return imp->GetAngleTag(t); }
    Dihedral GetDihedralTag(uid_ t) const { return imp->GetDihedralTag(t); }
    Angle GetAngleID(uid_ i) const { return imp->GetAngleID(i); }
    Dihedral GetDihedralID(uid_ i) const { return imp->GetDihedralID(i); }
    Atom GetAtom(size_ p) const { return imp->GetAtom(p); }
    Atom GetAtomTag(uid_ t) const { return imp->GetAtomTag(t); }
    Atom GetAtomID(uid_ i) const { return imp->GetAtomID(i); }
    Bond GetBond(size_ p) const { return imp->GetBond(p); }
    Bond GetBond(const Atom& a, const Atom& b) const { return imp->GetBond(a,b); }
    Bond GetBondTag(uid_ t) const { return imp->GetBondTag(t); }
    Bond GetBondID(uid_ i) const { return imp->GetBondID(i); }
    string_ GetFormula() { return imp->GetFormula(); }
    const graph::MolecularGraph& GetGraph() const { return imp->GetGraph(); }
    string_ GetName() const { return imp->GetName(); }
    int_ GetMolecularCharge() const { return imp->GetMolecularCharge(); }
    size_ NumAtoms() const { return imp->NumAtoms(); }
    size_ NumBonds() const { return imp->NumBonds(); }
    size_ NumAngles() { return imp->NumAngles(); }
    size_ NumDihedrals() { return imp->NumDihedrals(); }
    void SetName(string_ n) { imp->SetName(n); }
    void SetMolecularCharge(int_ q) { imp->SetMolecularCharge(q); }
    bool HasAtom(const Atom& a) const { return imp->HasAtom(a); }
    bool HasBond(const Bond& b) const { return imp->HasBond(b); }
    bool HasBond(const Atom& a, const Atom& b) const { return imp->HasBond(a,b); }
    bool HasAngle(const Angle& a) const { return imp->HasAngle(a); }
    bool HasAngle(const Atom& a, const Atom& b, const Atom& c) {
      return imp->HasAngle(a,b,c);
    }
    bool HasDihedral(const Dihedral& d) const { return imp->HasDihedral(d); }
    bool HasDihedral(const Atom& a, const Atom& b, const Atom& c,
                     const Atom& d) {
      return imp->HasDihedral(a,b,c,d);
    }
    Atom NewAtom() { return imp->NewAtom(); }
    Atom NewAtom(Element e) { return imp->NewAtom(e); }
    Atom NewAtom(string_ n) { return imp->NewAtom(n); }
    Atom NewAtom(string_ n, Element e) { return imp->NewAtom(n,e); }
    Bond NewBond(Atom a, Atom b) { return imp->NewBond(a,b); }
    bool RemoveAtom(Atom a) { return imp->RemoveAtom(a); }
    bool RemoveBond(Bond b) { return imp->RemoveBond(b); }
    bool RemoveBond(Atom a, Atom b) { return imp->RemoveBond(a,b); }
    size_ PerceiveAngles() { return imp->PerceiveAngles(); }
    size_ PerceiveDihedrals() { return imp->PerceiveDihedrals(); }
    void ReserveAtoms(size_ n) { return imp->ReserveAtoms(n); }
    void ReserveBonds(size_ n) { return imp->ReserveBonds(n); }
    void SetPropertyModified(MolProperty p ) { imp->SetPropertyModified(p); }
    Atoms GetAtoms() const { return imp->GetAtoms(); }
    Bonds GetBonds() const { return imp->GetBonds(); }
    Angles GetAngles() const { return imp->GetAngles(); }
    Dihedrals GetDihedrals() const { return imp->GetDihedrals(); }
    uid_ GetUniqueID() const { return imp->GetUniqueID(); }
    
    // Internals access
    string_ get_name() const { return imp->_name; }
    int_ get_q() const { return imp->_q; }
    IAtoms& get_atms() const { return imp->_atms; }
    IBonds& get_bnds() const { return imp->_bnds; }
    IAngles& get_angs() const { return imp->_angs; }
    IDihedrals& get_dhds() const { return imp->_dhds; }
    EmergeSet& get_emerge() const { return imp->_emerge; }
    graph::MolecularGraph get_g() const { return imp->_g; }
    string_ get_formula_cache() const { return imp->_formula_cache; }
  };
  
  struct MoleculeTestFixture {
    TestMolecule benzene, randmol, blankmol;
    std::vector<Atom> a_benzene, a_randmol;
    std::vector<Bond> b_benzene, b_randmol;
    std::vector<std::pair<size_, size_>> e_benzene, e_randmol;
    std::vector<stdx::triple<size_, size_, size_>> g_benzene;
    std::vector<stdx::quad<size_, size_, size_, size_>> d_benzene;
    std::random_device rd;
    std::mt19937 generator;
    
    using EmergeSet = TestMolecule::EmergeSet;
    using EmergeProp = TestMolecule::EmergeProp;
    
    MoleculeTestFixture() : benzene(true), randmol(true), blankmol(),
    generator(rd()) { }
    
    void BuildBenzene() {
      e_benzene.assign({{0,1},{0,5},{0,6},{1,7},{1,2},{2,3},{2,8},{3,4},
        {3,9},{4,5},{4,10},{5,11}});
      g_benzene.assign({{1,0,5},{1,0,6},{5,0,6},{7,1,0},{0,1,2},{7,1,2},
                        {1,2,3},{1,2,8},{3,2,8},{2,3,4},{2,3,9},{4,3,9},
                        {3,4,5},{3,4,10},{5,4,10},{0,5,4},{0,5,11},
                        {4,5,11}});
      d_benzene.assign({{7,1,0,5},{7,1,0,6},{2,1,0,5},{2,1,0,6},
                        {1,0,5,4},{1,0,5,11},{6,0,5,4},{6,0,5,11},
                        {0,1,2,3},{0,1,2,8},{7,1,2,3},{7,1,2,8},
                        {1,2,3,4},{1,2,3,9},{8,2,3,4},{8,2,3,9},
                        {2,3,4,5},{2,3,4,10},{9,3,4,5},{9,3,4,10},
                        {3,4,5,0},{3,4,5,11},{10,4,5,0},{10,4,5,11}});
      a_benzene.reserve(12);
      b_benzene.reserve(12);
      for (size_ i = 0; i < 12; ++i) {
        a_benzene.emplace_back(benzene.NewAtom());
        a_benzene.back()->SetTag(10 + i);
        if (i < 6) a_benzene.back()->SetElement("C");
        else if (i == 6) a_benzene.back()->SetElement("Cl");
        else if (i == 7) a_benzene.back()->SetElement("Br");
        else if (i == 8) a_benzene.back()->SetElement("I");
        else if (i == 9) a_benzene.back()->SetElement("F");
        else a_benzene.back()->SetElement("H");
      }
      for (auto& i : e_benzene) {
        b_benzene.emplace_back(benzene.NewBond(a_benzene[i.first],
                                                   a_benzene[i.second]));
        b_benzene.back()->SetTag(b_benzene.size());
      }
      b_benzene.front()->SwapSourceTarget();
      b_benzene.back()->SwapSourceTarget();
    }
    
    void BuildRandomMolecule() {
      std::uniform_int_distribution<size_> atom_count(15,35);
      std::uniform_int_distribution<size_> bond_count(25,55);
      std::uniform_int_distribution<size_> element(1, 50);
      size_ num_atoms = atom_count(generator);
      size_ num_bonds = bond_count(generator);
      a_randmol.reserve(num_atoms);
      b_randmol.reserve(num_bonds);
      randmol.ReserveAtoms(num_atoms);
      randmol.ReserveBonds(num_bonds);
      for (size_ i = 0; i < num_atoms; ++i) {
        a_randmol.emplace_back(randmol.NewAtom());
        a_randmol.back()->SetTag(i + num_atoms + 11);
        a_randmol.back()->SetElement(element(generator));
      }
      
      std::vector<std::pair<size_, size_>> possible_bonds;
      possible_bonds.reserve(num_atoms * (num_atoms + 1) / 2);
      for (size_ i = 0; i < num_atoms; ++i) {
        for (size_ j = i + 1; j < num_atoms; ++j )
          possible_bonds.emplace_back(i,j);
      }
      std::sample(possible_bonds.begin(), possible_bonds.end(),
                  std::back_inserter(e_randmol), num_bonds, generator);
      for (auto& i : e_randmol){
        b_randmol.emplace_back(randmol.NewBond(a_randmol[i.first],
                                                   a_randmol[i.second]));
        b_randmol.back()->SetTag(b_randmol.size() + num_bonds + 2);
      }
    }
  };
}

#endif /* INDIGOX_TEST_MOLECULE_TEST_HPP */
