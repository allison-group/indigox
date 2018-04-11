#include <indigox/python/pickle.hpp>

namespace py = pybind11;
using namespace indigox;

size_t __indigox_pickling_version = 0;

py::tuple PickleAtom(const Atom atom) {
  return py::make_tuple(atom->GetAromaticity(),
                        atom->GetElement()->GetAtomicNumber(),
                        atom->GetFormalCharge(),
                        atom->GetImplicitCount(),
                        atom->GetIndex(),
//                        atom->GetMolecule(),
                        atom->GetName(),
//                        atom->GetStereochemistry(),
//                        atom->GetVector(),
                        __indigox_pickling_version);
}

py::tuple PicklePeriodicTable(const PeriodicTable table) {
  return py::make_tuple(bool(table));
}

Atom __unpickle_atom_version_0(py::tuple& t) {
  Atom a = Atom(new IXAtom());
  a->SetAromaticity(t[0].cast<bool>());
  a->SetElement(t[1].cast<unsigned int>());
  a->SetFormalCharge(t[2].cast<int>());
  a->SetImplicitCount(t[3].cast<unsigned int>());
  a->SetIndex(t[4].cast<unsigned int>());
  a->SetName(t[5].cast<std::string>());
  return a;
}

Atom UnpickleAtom(py::tuple& t) {
  size_t pickled_version = t[t.size() - 1].cast<size_t>();
  switch (pickled_version) {
    case 0:
      return __unpickle_atom_version_0(t);
      
    default:
      return __unpickle_atom_version_0(t);
  }
}

PeriodicTable UnpicklePeriodicTable(py::tuple& t) {
  bool real = t[0].cast<bool>();
  return real ? IXPeriodicTable::GetInstance() : PeriodicTable();
}
