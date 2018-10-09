#ifndef INDIGOX_TEST_ATOM_TEST_HPP
#define INDIGOX_TEST_ATOM_TEST_HPP

#include <random>
#include <Eigen/Dense>

namespace indigox::test {
  struct TestAtom {
    // Typedefs
    using IAngles = IXAtom::AtomAngles;
    using IBonds = IXAtom::AtomBonds;
    using IDihedrals = IXAtom::AtomDihedrals;
    using Angles = std::pair<IXAtom::AtomAngleIter, IXAtom::AtomAngleIter>;
    using Bonds = std::pair<IXAtom::AtomBondIter, IXAtom::AtomBondIter>;
    using Dihedrals = std::pair<IXAtom::AtomDihedralIter, IXAtom::AtomDihedralIter>;
    
    indigox::Atom imp;
    
    // Private wrapping members
    TestAtom() = delete;
    TestAtom(Molecule m) : imp(new IXAtom(m)) { }
    void AddAngle(Angle a) { imp->AddAngle(a); }
    void AddBond(Bond b) { imp->AddBond(b); }
    void AddDihedral(Dihedral d) { imp->AddDihedral(d); }
    void Clear() { imp->Clear(); }
    void RemoveAngle(Angle a) { imp->RemoveAngle(a); }
    void RemoveBond(Bond b) { imp->RemoveBond(b); }
    void RemoveDihedral(Dihedral d) { imp->RemoveDihedral(d); }
    
    // Public wrapping members
    uint32_t GetUniqueID() const { return imp->GetUniqueID(); }
    Element GetElement() const { return imp->GetElement(); }
    int GetFormalCharge() const { return imp->GetFormalCharge(); }
    double GetPartialCharge() const { return imp->GetPartialCharge(); }
    uint32_t GetTag() const { return imp->GetTag(); }
    uint32_t GetImplicitCount() const { return imp->GetImplicitCount(); }
    uint32_t AddImplicitHydrogen() { return imp->AddImplicitHydrogen(); }
    uint32_t RemoveImplicitHydrogen() { return imp->RemoveImplicitHydrogen(); }
    Molecule GetMolecule() const { return imp->GetMolecule(); }
    std::string GetName() const { return imp->GetName(); }
    double GetX() const { return imp->GetX(); }
    double GetY() const { return imp->GetY(); }
    double GetZ() const { return imp->GetZ(); }
    const Eigen::Vector3d& GetVector() const { return imp->GetVector(); }
    std::string ToString() const { return imp->ToString(); }
    void SetElement(Element e) { imp->SetElement(e); }
    void SetElement(std::string e) { imp->SetElement(e); }
    void SetElement(uint32_t e) { imp->SetElement(e); }
    void SetFormalCharge(int q) { imp->SetFormalCharge(q); }
    void SetPartialCharge(double q) { imp->SetPartialCharge(q); }
    void SetImplicitCount(uint32_t h) { imp->SetImplicitCount(h); }
    void SetTag(uint32_t t) { imp->SetTag(t); }
    void SetName(std::string n) { imp->SetName(n); }
    void SetX(double x) { imp->SetX(x); }
    void SetY(double y) { imp->SetY(y); }
    void SetZ(double z) { imp->SetZ(z); }
    void SetPosition(double x, double y, double z) { imp->SetPosition(x,y,z); }
    void SetStereochemistry(AtomStereo s) { imp->SetStereochemistry(s); }
    void SetAromaticity(bool a) { imp->SetAromaticity(a); }
    AtomStereo GetStereochemistry() const { return imp->GetStereochemistry(); }
    bool GetAromaticity() const { return imp->GetAromaticity(); }
    Bonds GetBondIters() const { return imp->GetBondIters(); }
    Angles GetAngleIters() const { return imp->GetAngleIters(); }
    Dihedrals GetDihedralIters() const { return imp->GetDihedralIters(); }
    size_t NumBonds() const { return imp->NumBonds(); }
    size_t NumAngles() const { return imp->NumAngles(); }
    size_t NumDihedrals() const { return imp->NumDihedrals(); }
    size_t GetIndex() const { return imp->GetIndex(); }
    FFAtom GetType() const { return imp->GetType(); }
    void SetType(FFAtom t) { imp->SetType(t); }

    // Internals access
    _Molecule get_mol() const { return imp->_mol; }
    _Element get_elem() const { return imp->_elem; }
    int get_fc() const { return imp->_fc; }
    uint32_t get_tag() const { return imp->_tag; }
    uint32_t get_implicitH() const { return imp->_implicitH; }
    std::string get_name() const { return imp->_name; }
    const Eigen::Vector3d& get_pos() const { return imp->_pos; }
    double get_partial() const { return imp->_partial; }
    AtomStereo get_stereo() const { return imp->_stereo; }
    bool get_aromatic() const { return imp->_aromatic; }
    const IBonds& get_bnds() const { return imp->_bnds; }
    const IAngles& get_angs() const { return imp->_angs; }
    const IDihedrals& get_dhds() const { return imp->_dhds; }
    FFAtom get_type() const { return imp->_type; }
  };
  
  inline TestAtom CreateGenericTestAtom() {
    return TestAtom(Molecule());
  }
  
  struct AtomTestFixture {
    Molecule mol = CreateMolecule();
    TestAtom atm;
    FFAtom fftype;
    AtomTestFixture();
  };
}

#endif /* INDIGOX_TEST_ATOM_TEST_HPP */
