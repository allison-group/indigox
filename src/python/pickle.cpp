#include <indigox/python/pickle.hpp>

namespace py = pybind11;
using namespace indigox;

size_t __indigox_pickling_version = 0;

py::tuple PickleAtom(const Atom a) {
  return py::make_tuple(a->GetAromaticity(),
                        a->GetElement()->GetAtomicNumber(),
                        a->GetFormalCharge(),
                        a->GetImplicitCount(),
                        a->GetIndex(),
//                        a->GetMolecule(),
                        a->GetName(),
//                        a->GetStereochemistry(),
//                        a->GetVector(),
                        __indigox_pickling_version);
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
