#ifndef INDIGOX_TEST_PERIODICTABLE_TEST_HPP
#define INDIGOX_TEST_PERIODICTABLE_TEST_HPP

namespace indigox::test {
  struct TestPeriodicTable {
    // Typedefs
    using ZType = const indigox::IXPeriodicTable::atomic_number_map_type &;
    using XType = const indigox::IXPeriodicTable::atomic_symbol_map_type &;
    using EType = indigox::IXPeriodicTable::element_type;
    
    indigox::PeriodicTable imp;
    
    // Private wrapping functions
    TestPeriodicTable() : imp(new IXPeriodicTable()) { }
    void GeneratePeriodicTable() { imp->GeneratePeriodicTable(); }
    
    // Public wrapping functions
    EType GetElement(const uchar_ z) const { return imp->GetElement(z); }
    EType GetElement(const string_ n) const { return imp->GetElement(n); }
    EType operator[](const uchar_ z) const { return (*imp)[z]; }
    EType operator[](const string_ n) const { return (*imp)[n]; }
    EType GetUndefined() const { return imp->GetUndefined(); }
    size_ NumElements() const { return imp->NumElements(); }
    string_ ToString() const { return imp->ToString(); }
    
    // Internals access
    EType get_null() const { return imp->_null; }
    ZType get_z_to() const { return imp->_z_to; }
    XType get_name_to() const { return imp->_name_to; }
    
  };
  
  struct TestElement {
    indigox::Element imp;
    using u = uchar_;
    using s = string_;
    using f = float_;
    // Private wrapping functions
    TestElement() = delete;
    TestElement(u Z, s N, s S, f M, u G, u P, u V, u O, u H, f R, f C, f W, f X)
    : imp(new IXElement(Z, N, S, M, G, P, V, O, H, R, C, W, X)) { }
    
    // Public wrapping functions
    float_ GetAtomicMass() const { return imp->GetAtomicMass(); }
    size_ GetAtomicNumber() const { return imp->GetAtomicNumber(); }
    float_ GetAtomicRadius() const { return imp->GetAtomicRadius(); }
    float_ GetCovalentRadius() const { return imp->GetCovalentRadius(); }
    float_ GetVanDerWaalsRadius() const { return imp->GetVanDerWaalsRadius(); }
    string_ GetName() const { return imp->GetName(); }
    string_ GetSymbol() const { return imp->GetSymbol(); }
    uchar_ GetGroup() const { return imp->GetGroup(); }
    uchar_ GetPeriod() const { return imp->GetPeriod(); }
    uchar_ GetValenceElectronCount() const { return imp->GetValenceElectronCount(); }
    uchar_ GetOctet() const { return imp->GetOctet(); }
    uchar_ GetHypervalentOctet() const { return imp->GetHypervalentOctet(); }
    float_ GetElectronegativity() const { return imp->GetElectronegativity(); }
    string_ ToString() const { return imp->ToString(); }
    
    // Internals access
    string_ get_nme() const { return imp->_nme; }
    string_ get_sym() const { return imp->_sym; }
    uchar_ get_grp() const { return imp->_grp; }
    uchar_ get_prd() const { return imp->_prd; }
    uchar_ get_Z() const { return imp->_Z; }
    uchar_ get_val() const { return imp->_val; }
    uchar_ get_oct() const { return imp->_oct; }
    uchar_ get_hyp() const { return imp->_hyp; }
    float_ get_mass() const { return imp->_mass; }
    float_ get_rad() const { return imp->_rad; }
    float_ get_cov() const { return imp->_cov; }
    float_ get_vdw() const { return imp->_vdw; }
    float_ get_chi() const { return imp->_chi; }
  };
  
  inline TestElement CreateGenericTestElement() {
    return TestElement(23, "Testium", "Tm", 12.34, 5, 6, 7, 8, 9, 10.11, 12.13,
                       14.1789, 15.8763);
  }
}

#endif // INDIGOX_TEST_PERIODICTABLE_TEST_HPP
