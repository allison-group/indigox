#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/json.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace indigox {
  // ===========================================================================
  // == Element Data implementation ============================================
  // ===========================================================================

  struct Element::ElementImpl {
    std::string name, symbol;
    int32_t group, period, atomic_number, valence, octet, hyper_octet;
    double mass, radius, vanderwaals, chi;
    std::array<double, 3> covalent;
    ElementImpl() : name("Undefined"), symbol("XX"), group(-1), period(-1), atomic_number(-1), valence(-1), octet(-1), hyper_octet(-1), mass(-1), radius(-1), vanderwaals(-1), chi(-1) {
      covalent.fill(-1);
    }
    ElementImpl(std::string nme, std::string sym, int32_t grp, int32_t prd,
                int32_t num, int32_t val, int32_t oct, int32_t hyp, double mas,
                double rad, std::array<double, 3>& cov, double vdw, double eln)
        : name(nme), symbol(sym), group(grp), period(prd), atomic_number(num),
          valence(val), octet(oct), hyper_octet(hyp), mass(mas), radius(rad),
           vanderwaals(vdw), chi(eln), covalent(cov) {}

    static std::shared_ptr<ElementImpl> GetNullState() {
      static std::shared_ptr<ElementImpl> state;
      if (!state)
        state = std::make_shared<ElementImpl>();
      return state;
    }
  };

  // ===========================================================================
  // == Element Construction and Assignment ====================================
  // ===========================================================================

  Element::Element() : m_data(ElementImpl::GetNullState()) {}
  Element::Element(const Element &e2) : m_data(e2.m_data) {}
  Element::Element(Element &&e2) : m_data(std::move(e2.m_data)) {}
  Element &Element::operator=(const Element &e2) {
    if (&e2 != this) m_data = e2.m_data;
    return *this;
  }
  Element &Element::operator=(Element &&e2) {
    m_data = std::move(e2.m_data);
    return *this;
  }
  Element::Element(uint8_t Z, std::string name, std::string sym, double mass,
                   uint8_t grp, uint8_t prd, uint8_t val, uint8_t oct,
                   uint8_t hyp, double rad, std::array<double, 3>& cov, double vdw, double chi)
      : m_data(std::make_shared<ElementImpl>(
            name, sym, grp, prd, Z, val, oct, hyp, mass, rad, cov, vdw, chi)) {}

  // ===========================================================================
  // == Element Data Retrevial =================================================
  // ===========================================================================

  double Element::GetAtomicMass() const { return m_data->mass; }
  int32_t Element::GetAtomicNumber() const { return m_data->atomic_number; }
  double Element::GetAtomicRadius() const { return m_data->radius; }
  double Element::GetCovalentRadius() const { return m_data->covalent[0]; }
  double Element::GetDoubleCovalentRadius() const { return m_data->covalent[1]; }
  double Element::GetTripleCovalentRadius() const { return m_data->covalent[2]; }
  double Element::GetVanDerWaalsRadius() const {
    return m_data->vanderwaals;
  }
  std::string Element::GetName() const { return m_data->name; }
  std::string Element::GetSymbol() const { return m_data->symbol; }
  int32_t Element::GetGroup() const { return m_data->group; }
  int32_t Element::GetPeriod() const { return m_data->period; }
  int32_t Element::GetValenceElectronCount() const {
    return m_data->valence;
  }
  int32_t Element::GetOctet() const { return m_data->octet; }
  int32_t Element::GetHypervalentOctet() const {
    return m_data->hyper_octet;
  }
  double Element::GetElectronegativity() const { return m_data->chi; }

  // ===========================================================================
  // == Element Operators ======================================================
  // ===========================================================================

  Element::operator bool() { return m_data->atomic_number != -1; }
  std::ostream &operator<<(std::ostream &os, const Element &element) {
    return (os << element.GetSymbol());
  }
  bool Element::operator==(int32_t Z) const { return Z == GetAtomicNumber(); }
  bool Element::operator==(const std::string &name) const {
    return (name.size() <= 2 ? name == GetSymbol()
                             : utils::ToUpperFirst(name) == GetName());
  }
  bool Element::operator==(const Element &element) const {
    return element.m_data == m_data;
  }
  bool Element::operator!=(int32_t Z) const { return !(*this == Z); }
  bool Element::operator!=(const std::string &Z) const { return !(*this == Z); }
  bool Element::operator!=(const Element &Z) const { return !(*this == Z); }
  bool Element::operator<(const Element &element) const {
    return GetAtomicNumber() < element.GetAtomicNumber();
  }
  bool Element::operator>(const Element &element) const {
    return GetAtomicNumber() > element.GetAtomicNumber();
  }
  bool Element::operator<=(const Element &element) const {
    return GetAtomicNumber() <= element.GetAtomicNumber();
  }
  bool Element::operator>=(const Element &element) const {
    return GetAtomicNumber() >= element.GetAtomicNumber();
  }

  // ===========================================================================
  // == PeriodicTable Construction/Generation ==================================
  // ===========================================================================

  PeriodicTable::PeriodicTable() {}

  void PeriodicTable::GeneratePeriodicTable() {
    using json = nlohmann::json;

    std::string pt_path = std::string(IX_DATA_DIRECTORY);
    if (pt_path.back() != '/') pt_path.append("/");
    pt_path.append("periodictable.json");

    std::ifstream pt_file(pt_path);
    json pt_dat;
    pt_file >> pt_dat;

    // Check that the json data is 119 items with the first being
    if (pt_dat.size() != 119) throw std::runtime_error("Expected 119 items in periodic table file.");
    
    _elems[0] = Element();
    
    int32_t count = 1;
    for (json::iterator begin = pt_dat.begin() + 1; begin != pt_dat.end(); ++begin, ++count) {
      std::vector<std::string> required_keys = {"Z", "atomicradius", "covalentradius_1", "covalentradius_2", "covalentradius_3", "electronegativity", "group", "mass", "name", "octet_hypervalent", "octet_normal", "period", "symbol", "valence", "vanderwaalsradius"}, optional_keys;
      
      std::string name, symbol, error_str;
      int32_t group, period, atomic_number, valence, octet, hyper_octet;
      double mass, radius, vanderwaals, chi;
      std::array<double, 3> covalent;
      error_str = utils::JSONKeyChecker(required_keys, optional_keys, *begin);
      if (!error_str.empty()) throw std::runtime_error(error_str);
      json data = *begin;
      
      name = data["name"];
      symbol = data["symbol"];
      group = data["group"];
      period = data["period"];
      atomic_number = data["Z"];
      valence = data["valence"];
      octet = data["octet_normal"];
      hyper_octet = data["octet_hypervalent"];
      mass = data["mass"];
      radius = data["atomicradius"];
      vanderwaals = data["vanderwaalsradius"];
      chi = data["electronegativity"];
      covalent[0] = data["covalentradius_1"];
      covalent[1] = data["covalentradius_2"];
      covalent[2] = data["covalentradius_3"];
      
      if (atomic_number != count) throw std::runtime_error("Unexpected element number");
      _elems[atomic_number] = Element(atomic_number, name, symbol, mass, group, period, valence, octet, hyper_octet, radius, covalent, vanderwaals, chi);
    }
    
    for (Element &e : _elems) {
      _name_to_idx.emplace(e.GetSymbol(), e.GetAtomicNumber());
      _name_to_idx.emplace(e.GetName(), e.GetAtomicNumber());
    }
  }

  // ===========================================================================
  // == PeriodicTable Data Retrevial ===========================================
  // ===========================================================================

  Element PeriodicTable::GetElement(const int32_t z) const {
    return (z > 0 && z < 119) ? _elems.at(z) : _elems.at(0);
  }

  Element PeriodicTable::GetElement(std::string name) const {
    if (name.size() > 2) name = utils::ToUpperFirst(name);
    auto pos = _name_to_idx.find(name);
    return (pos == _name_to_idx.end()) ? _elems.at(0) : _elems.at(pos->second);
  }

  // ===========================================================================
  // == PeriodicTable Operators ================================================
  // ===========================================================================

  using erow = std::pair<const Element &, int>;
  std::ostream &operator<<(std::ostream &ss, erow er) {
    switch (er.second) {
    case 0: ss << "--- "; break;
    case 1: ss << std::setw(3) << er.first.GetAtomicNumber() << '|'; break;
    case 2: ss << std::setw(3) << er.first.GetSymbol() << '|'; break;
    case 3: ss << "    "; break;
    case 4:
    default: ss << "   |"; break;
    }
    return ss;
  }

  std::string PeriodicTable_ToString() {
    const PeriodicTable &PT = GetPeriodicTable();
    std::stringstream ss;
    size_t row_count = 0, restart = 0;
    std::vector<int> elems = {
        1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   2,   -1,  3,   4,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   5,   6,   7,   8,   9,   10,  -1,  11,  12,  0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   13,  14,  15,  16,  17,  18,
        -1,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
        32,  33,  34,  35,  36,  -1,  37,  38,  39,  40,  41,  42,  43,  44,
        45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  -1,  55,  56,  57,
        72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,
        86,  -1,  87,  88,  89,  104, 105, 106, 107, 108, 109, 110, 111, 112,
        113, 114, 115, 116, 117, 118, -1,  0,   0,   0,   58,  59,  60,  61,
        62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  0,   -1,  0,   0,
        0,   90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101, 102,
        103, 0,   -1};
    if (!PT.NumElements()) return "Empty PeriodicTable";
    ss << ' ';
    const Element &_null = PT.GetUndefined();
    for (size_t i = 0; i < elems.size();) {
      if (elems[i] == -1) {
        ss << "\n";
        if ((row_count == 0 || row_count == 1) && elems[i - 18] != 0)
          ss << '|';
        else
          ss << ' ';
        ++row_count;
        if (row_count < 3)
          i = restart;
        else {
          ++i;
          restart = i;
          row_count = 0;
        }
      } else if (elems[i] == 0 && i > 19 && elems[i - 19] != 0 &&
                 row_count == 0) {
        ss << erow(_null, 0);
        ++i;
      } else if (elems[i] == 0 && elems[i + 1] == 0) {
        ss << erow(_null, 3);
        ++i;
      } else if (elems[i] == 0 && elems[i + 1] > 0 && row_count != 0) {
        ss << erow(_null, 4);
        ++i;
      } else if (elems[i] == 0 && elems[i + 1] > 0 && row_count == 0) {
        ss << erow(_null, 3);
        ++i;
      } else if (elems[i] == 0 && elems[i + 1] == -1) {
        ss << erow(_null, 3);
        ++i;
      } else {
        ss << erow(PT[elems[i]], row_count);
        ++i;
      }
    }
    // Final line at bottom of table
    for (size_t i = 0; i < 17; ++i) {
      if (i < 3)
        ss << erow(_null, 3);
      else
        ss << erow(_null, 0);
    }
    return ss.str();
  }

  std::ostream &operator<<(std::ostream &os, const PeriodicTable &) {
    return (os << PeriodicTable_ToString());
  }

  // ===========================================================================
  // == PeriodicTable External access ==========================================
  // ===========================================================================

  const PeriodicTable &GetPeriodicTable() {
    using pPeriodicTable = std::unique_ptr<PeriodicTable>;
    static pPeriodicTable instance = pPeriodicTable();
    if (!instance) {
      instance.reset(new PeriodicTable());
      instance->GeneratePeriodicTable();
    }
    return *instance;
  }

  // ===========================================================================
  // == PeriodicTable and Element Testing ======================================
  // ===========================================================================

  /*  test_case("IXElement construction") {
    test::TestElement e = test::CreateGenericTestElement();
    check_eq(e.get_Z(), 23);
    check_eq(e.get_nme(), "Testium");
    check_eq(e.get_sym(), "Tm");
    check_eq(e.get_mass(), approximately(12.34));
    check_eq(e.get_grp(), 5);
    check_eq(e.get_prd(), 6);
    check_eq(e.get_val(), 7);
    check_eq(e.get_oct(), 8);
    check_eq(e.get_hyp(), 9);
    check_eq(e.get_rad(), approximately(10.11));
    check_eq(e.get_cov(), approximately(12.13));
    check_eq(e.get_vdw(), approximately(14.1789));
    check_eq(e.get_chi(), approximately(15.8763));
    }
    */

  /*  test_case("IXElement printing methods") {
    test::TestElement e = test::CreateGenericTestElement();
    check_eq(e.ToString(), "Testium(23, Tm)");
    std::stringstream ss;
    ss.str(""); ss << e.imp;
    check_eq(ss.str(), "Element(Testium)");
    ss.str(""); ss << *e.imp;
    check_eq(ss.str(), "Element(Testium)");
    ss.str(""); ss << Element();
    check_eq(ss.str(), "");
    }

    test_case("IXElement comparison methods") {
    Element e1 = test::CreateGenericTestElement().imp;
    Element e2 = test::TestElement(6, "Car", "C", 0.0, 5, 6, 7, 8, 9, 0.1,
    2.3, 4.9, 5.3).imp;

    // Check both shared_ptr comparisons
    check_ne(e1, e2);
    check_eq(e1, e1);
    // Check one direct, one shared_ptr comparison
    check_eq(e1, *e1);
    check_eq(*e1, e1);
    check_ne(e2, *e1);
    check_ne(*e2, e1);
    // Check both direct comparisons
    check_ne(*e1, *e2);
    check_eq(*e1, e1);
    // Check shared_ptr comparisions with integers (atomic number)
    check_eq(e1, 23);
    check_eq(23, e1);
    check_ne(23, e2);
    check_ne(e2, 23);
    // Check direct comparisions with integers (atomic number)
    check_eq(*e2, 6);
    check_eq(6, *e2);
    check_ne(6, *e1);
    check_ne(*e1, 6);
    // Check shared_ptr comparisons with strings (name)
    check_eq(e1, "Testium");
    check_eq(e1, "tesTIUM"); // name compare should be case insensitive
    check_eq("Testium", e1);
    check_eq("CAR", e2);     // name compare should be case insensitive
    check_ne("FAIL", e1);
    check_ne(e1, "FAIL");
    // Check shared_ptr comparisons with strings (symbol)
    check_eq(e1, "Tm");
    check_ne(e2, "Tm");
    check_ne(e1, "tm");  // symbol compare shoud be case sensitive
    check_eq("C", e2);
    check_ne("c", e2);   // symbol compare should be case sensitive
    // Check direct comparisons with strings (name)
    check_eq(*e1, "Testium");
    check_eq(*e1, "tesTIUM"); // name compare should be case insensitive
    check_eq("Testium", *e1);
    check_eq("CAR", *e2);
    check_ne("FAIL", *e1);
    check_ne(*e1, "FAIL");
    // Check direct comparisons with strings (symbol)
    check_eq(*e1, "Tm");
    check_ne(*e2, "Tm");
    check_ne(*e1, "tm");  // symbol compare shoud be case sensitive
    check_eq("C", *e2);
    check_ne("c", *e2);
    }

    test_case("IXElement property getting") {
    test::TestElement e = test::CreateGenericTestElement();
    check_eq(e.GetAtomicNumber(), 23);
    check_eq(e.GetName(), "Testium");
    check_eq(e.GetSymbol(), "Tm");
    check_eq(e.GetAtomicMass(), approximately(12.34));
    check_eq(e.GetGroup(), 5);
    check_eq(e.GetPeriod(), 6);
    check_eq(e.GetValenceElectronCount(), 7);
    check_eq(e.GetOctet(), 8);
    check_eq(e.GetHypervalentOctet(), 9);
    check_eq(e.GetAtomicRadius(), approximately(10.11));
    check_eq(e.GetCovalentRadius(), approximately(12.13));
    check_eq(e.GetVanDerWaalsRadius(), approximately(14.1789));
    check_eq(e.GetElectronegativity(), approximately(15.8763));
    }
    */

  /*  test_case("IXPeriodicTable printing methods") {
      std::stringstream os;
      os << " --- --- \n"; os << "|  1| |  2|\n"; os << "|  H| | He|\n"; os << "
    --- ---                                         --- --- --- --- --- --- \n";
      os << "|  3|  4|                                       |  5|  6|  7|  8|
    9| 10|\n"; os << "| Li| Be|                                       |  B|  C|
    N|  O|  F| Ne|\n"; os << " --- --- --- --- --- --- --- --- \n"; os << "| 11|
    12|                                       | 13| 14| 15| 16| 17| 18|\n"; os
    << "| Na| Mg|                                       | Al| Si|  P|  S| Cl|
    Ar|\n"; os << " --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    --- --- --- \n"; os << "| 19| 20| 21| 22| 23| 24| 25| 26| 27| 28| 29| 30|
    31| 32| 33| 34| 35| 36|\n"; os << "|  K| Ca| Sc| Ti|  V| Cr| Mn| Fe| Co| Ni|
    Cu| Zn| Ga| Ge| As| Se| Br| Kr|\n"; os << " --- --- --- --- --- --- --- ---
    --- --- --- --- --- --- --- --- --- --- \n"; os << "| 37| 38| 39| 40| 41|
    42| 43| 44| 45| 46| 47| 48| 49| 50| 51| 52| 53| 54|\n"; os << "| Rb| Sr|  Y|
    Zr| Nb| Mo| Tc| Ru| Rh| Pd| Ag| Cd| In| Sn| Sb| Te|  I| Xe|\n"; os << " ---
    --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n"; os
    << "| 55| 56| 57| 72| 73| 74| 75| 76| 77| 78| 79| 80| 81| 82| 83| 84| 85|
    86|\n"; os << "| Cs| Ba| La| Hf| Ta|  W| Re| Os| Ir| Pt| Au| Hg| Tl| Pb| Bi|
    Po| At| Rn|\n"; os << " --- --- --- --- --- --- --- --- --- --- --- --- ---
    --- --- --- --- --- \n"; os << "| 87| 88|
    89|104|105|106|107|108|109|110|111|112|113|114|115|116|117|118|\n"; os << "|
    Fr| Ra| Ac| Db| Jl| Rf| Bh| Hn| Mt| Ds| Rg| Cn| Nh| Fl| Mc| Lv| Ts| Og|\n";
      os << " --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    --- --- \n"; os << "            | 58| 59| 60| 61| 62| 63| 64| 65| 66| 67|
    68| 69| 70| 71|    \n"; os << "            | Ce| Pr| Nd| Pm| Sm| Eu| Gd| Tb|
    Dy| Ho| Er| Tm| Yb| Lu|    \n"; os << "             --- --- --- --- --- ---
    --- --- --- --- --- --- --- ---     \n"; os << "            | 90| 91| 92|
    93| 94| 95| 96| 97| 98| 99|100|101|102|103|    \n"; os << "            | Th|
    Pa|  U| Np| Pu| Am| Cm| Bk| Cf| Es| Fm| Md| No| Lr|    \n"; os << " --- ---
    --- --- --- --- --- --- --- --- --- --- --- --- "; std::string
    expected_to_string = os.str();

      test::TestPeriodicTable pt;
      os.str(""); os << PeriodicTable();
      check_eq(os.str(), "");
      os.str(""); os << pt.imp;
      check_eq(os.str(), "PeriodicTable(0 elements)");
      os.str(""); os << *pt.imp;
      check_eq(os.str(), "PeriodicTable(0 elements)");
      check_eq("Empty PeriodicTable", pt.ToString());

      pt.GeneratePeriodicTable();
      os.str(""); os << pt.imp;
      check_eq(os.str(), "PeriodicTable(118 elements)");
      os.str(""); os << *pt.imp;
      check_eq(os.str(), "PeriodicTable(118 elements)");
      check_eq(expected_to_string, pt.ToString());
    }
    */

  /*  test_case("IXPeriodicTable getting elements") {
      // We build our test elements to compare the gets with
      test::TestPeriodicTable PT; PT.GeneratePeriodicTable();
      Element C = test::TestElement(6, "Carbon" , "C", 12.011000, 14, 2, 4, 8,
    8, 0.77, 0.77, 1.85, 2.55).imp; Element K = test::TestElement(19,
    "Potassium" , "K", 39.098300, 1, 4, 1, 2, 2, 2.27, 2.03, 2.31, 0.82).imp;
      Element null_element = PT.get_null();

      check_eq(null_element, PT.GetUndefined());
      // Correct retrival
      check_eq(C, PT.GetElement(6));
      check_eq(K, PT.GetElement(19));
      check_eq(C, PT.GetElement("C"));
      check_eq(K, PT.GetElement("K"));
      check_eq(C, PT.GetElement("Carbon"));
      check_eq(K, PT.GetElement("Potassium"));
      // Symbol get case sensitive, name case insensitive
      check_ne(C, PT.GetElement("c"));
      check_eq(K, PT.GetElement("potASSIUM"));
      // Bad requests return null atom
      check_eq(null_element, PT.GetElement("k"));
      check_eq(null_element, PT.GetElement("NotAnElement"));
      check_eq(null_element, PT.GetElement(0));
      check_eq(null_element, PT.GetElement(255));
      // operators work
      check_eq(C, PT["C"]);
      check_eq(K, PT[19]);
      check_eq(C, PT["Carbon"]);
      check_eq(null_element, PT[255]);
      check_eq(null_element, PT["NotAnElement"]);

    }
   */

  /*  test_case("IXPeriodicTable generation") {
      test::TestPeriodicTable PT; PT.GeneratePeriodicTable();
      test::TestPeriodicTable::ZType& z_to_element = PT.get_z_to();
      test::TestPeriodicTable::XType& name_to_element = PT.get_name_to();

      std::string name;
      uint8_t z;
      Element e;
      subcase("Checking atomic number table") {
        check_eq(z_to_element.size(), INDIGOX_NUM_ELEMENTS);
        DOCTEST_INFO("Current atomic number:");
        for (auto& ze : z_to_element) {
          std::tie(z, e) = ze;
          DOCTEST_CAPTURE(z);
          check_eq(e, name_to_element.at(e->GetSymbol()));
          check_eq(e, name_to_element.at(e->GetName()));
        }
      }

      subcase("Checking atomic name/symbol table") {
        check_eq(name_to_element.size(), INDIGOX_NUM_ELEMENTS * 2);
        DOCTEST_INFO("Current atomic symbol/name:");
        for (auto& ne : name_to_element) {
          std::tie(name, e) = ne;
          DOCTEST_CAPTURE(name);
          check_eq(e, z_to_element.at(e->GetAtomicNumber()));
        }
      }

      subcase("Checking null atom") {
        IXElement null_ = *PT.get_null();
        check_eq(null_.GetName(), "Undefined");
        check_eq(null_.GetSymbol(), "XX");
        check_eq(null_.GetAtomicNumber(), 0);
        check_eq(PT.NumElements(), INDIGOX_NUM_ELEMENTS);
      }

   }
   */

  /*  test_case("IXPeriodicTable GetPeriodicTable helper method") {
      PeriodicTable PT = GetPeriodicTable();
      check_eq(PT->NumElements(), INDIGOX_NUM_ELEMENTS);
      PeriodicTable PT2 = GetPeriodicTable();
      check_eq(PT2->NumElements(), INDIGOX_NUM_ELEMENTS);
      check_eq(PT.get(), PT2.get());
    }
   */

} // namespace indigox
