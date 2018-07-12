#include <cstdint>
#include <set>
#include <sstream>
#include <string>

#include <boost/algorithm/string/join.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <indigox/python/interface.hpp>
#include <indigox/classes/periodictable.hpp>

namespace py = pybind11;

void GeneratePyElement(py::module& m) {
  using namespace indigox;
  py::class_<IXElement, Element>(m, "Element")
  // No constructor
  .def("__eq__", py::overload_cast<Element, uint8_t>(&operator==))
  .def("__eq__", py::overload_cast<uint8_t, Element>(&operator==))
  .def("__eq__", py::overload_cast<Element, std::string>(&operator==))
  .def("__eq__", py::overload_cast<std::string, Element>(&operator==))
  .def("__eq__", py::overload_cast<Element, Element>(&operator==))
  .def("__ne__", py::overload_cast<Element, uint8_t>(&operator!=))
  .def("__ne__", py::overload_cast<uint8_t, Element>(&operator!=))
  .def("__ne__", py::overload_cast<Element, std::string>(&operator!=))
  .def("__ne__", py::overload_cast<std::string, Element>(&operator!=))
  .def("__ne__", py::overload_cast<Element, Element>(&operator!=))
  .def("__repr__", [](Element e) {
    std::stringstream ss; ss << e; return ss.str();
  }) // Data print for developers
  .def("GetAtomicMass", &IXElement::GetAtomicMass)
  .def("GetAtomicNumber", &IXElement::GetAtomicNumber)
  .def("GetAtomicRadius", &IXElement::GetAtomicRadius)
  .def("GetCovalentRadius", &IXElement::GetCovalentRadius)
  .def("GetElectronegativity", &IXElement::GetElectronegativity)
  .def("GetGroup", &IXElement::GetGroup)
  .def("GetHypervalentOctet", &IXElement::GetHypervalentOctet)
  .def("GetName", &IXElement::GetName)
  .def("GetOctet", &IXElement::GetOctet)
  .def("GetPeriod", &IXElement::GetPeriod)
  .def("GetSymbol", &IXElement::GetSymbol)
  .def("GetVanDerWaalsRadius", &IXElement::GetVanDerWaalsRadius)
  .def("GetValenceElectronCount", &IXElement::GetValenceElectronCount)
  .def("ToString", &IXElement::ToString)
  // No pickling support
  ;
  
  using ElementSet = std::set<Element>;
  py::class_<ElementSet>(m, "ElementSet")
  .def(py::init<>())
  .def("__len__", [](const ElementSet& s) { return s.size(); })
  .def("__iter__", [](ElementSet& s) {
    return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0, 1>())
  .def("__contains__", [](const ElementSet& s, const Element& e) {
    return bool(s.count(e)); })
  .def("__repr__", [](const ElementSet& s) {
    std::vector<string_> strs;
    strs.reserve(s.size());
    std::stringstream ss;
    for (Element e : s) {
      ss << e;
      strs.push_back(ss.str());
      ss.str("");
    }
    ss << "{" << boost::join(strs, ", ") << "}";
    return ss.str();
  })
  .def("clear", &ElementSet::clear)
  .def("insert", (std::pair<ElementSet::iterator, bool> (ElementSet::*)(const Element&)) &ElementSet::insert)
  .def("erase", (size_t (ElementSet::*)(const Element&)) &ElementSet::erase)
  ;
  
}

