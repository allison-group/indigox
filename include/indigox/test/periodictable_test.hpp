#ifndef INDIGOX_TEST_PERIODICTABLE_TEST_HPP
#define INDIGOX_TEST_PERIODICTABLE_TEST_HPP

#include <cstdint>
#include <string>

namespace indigox::test {
  struct TestPeriodicTable {
//    // Typedefs
//    using ZType = const indigox::IXPeriodicTable::atomic_number_map_type &;
//    using XType = const indigox::IXPeriodicTable::atomic_symbol_map_type &;
//    using EType = indigox::IXPeriodicTable::element_type;
//
//    indigox::PeriodicTable imp;
//
//    // Private wrapping functions
//    TestPeriodicTable() : imp(new IXPeriodicTable()) { }
//    void GeneratePeriodicTable() { imp->GeneratePeriodicTable(); }
//
//    // Public wrapping functions
//    EType GetElement(const uint8_t z) const { return imp->GetElement(z); }
//    EType GetElement(const std::string n) const { return imp->GetElement(n); }
//    EType operator[](const uint8_t z) const { return (*imp)[z]; }
//    EType operator[](const std::string n) const { return (*imp)[n]; }
//    EType GetUndefined() const { return imp->GetUndefined(); }
//    size_t NumElements() const { return imp->NumElements(); }
//    std::string ToString() const { return imp->ToString(); }
//
//    // Internals access
//    EType get_null() const { return imp->_null; }
//    ZType get_z_to() const { return imp->_z_to; }
//    XType get_name_to() const { return imp->_name_to; }
//
  };

  struct TestElement {
//    indigox::Element imp;
//    using u = uint8_t;
//    using s = std::string;
//    using f = double;
//    // Private wrapping functions
//    TestElement() = delete;
//    TestElement(u Z, s N, s S, f M, u G, u P, u V, u O, u H, f R, f C, f W, f X)
//    : imp(new IXElement(Z, N, S, M, G, P, V, O, H, R, C, W, X)) { }
//
//    // Public wrapping functions
//    double GetAtomicMass() const { return imp->GetAtomicMass(); }
//    size_t GetAtomicNumber() const { return imp->GetAtomicNumber(); }
//    double GetAtomicRadius() const { return imp->GetAtomicRadius(); }
//    double GetCovalentRadius() const { return imp->GetCovalentRadius(); }
//    double GetVanDerWaalsRadius() const { return imp->GetVanDerWaalsRadius(); }
//    std::string GetName() const { return imp->GetName(); }
//    std::string GetSymbol() const { return imp->GetSymbol(); }
//    uint8_t GetGroup() const { return imp->GetGroup(); }
//    uint8_t GetPeriod() const { return imp->GetPeriod(); }
//    uint8_t GetValenceElectronCount() const { return imp->GetValenceElectronCount(); }
//    uint8_t GetOctet() const { return imp->GetOctet(); }
//    uint8_t GetHypervalentOctet() const { return imp->GetHypervalentOctet(); }
//    double GetElectronegativity() const { return imp->GetElectronegativity(); }
//    std::string ToString() const { return imp->ToString(); }
//
//    // Internals access
//    std::string get_nme() const { return imp->_nme; }
//    std::string get_sym() const { return imp->_sym; }
//    uint8_t get_grp() const { return imp->_grp; }
//    uint8_t get_prd() const { return imp->_prd; }
//    uint8_t get_Z() const { return imp->_Z; }
//    uint8_t get_val() const { return imp->_val; }
//    uint8_t get_oct() const { return imp->_oct; }
//    uint8_t get_hyp() const { return imp->_hyp; }
//    double get_mass() const { return imp->_mass; }
//    double get_rad() const { return imp->_rad; }
//    double get_cov() const { return imp->_cov; }
//    double get_vdw() const { return imp->_vdw; }
//    double get_chi() const { return imp->_chi; }
  };
//
//  inline TestElement CreateGenericTestElement() {
//    return TestElement(23, "Testium", "Tm", 12.34, 5, 6, 7, 8, 9, 10.11, 12.13,
//                       14.1789, 15.8763);
//  }
}

#endif // INDIGOX_TEST_PERIODICTABLE_TEST_HPP
