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
                        atom->GetStereochemistry(),
//                        atom->GetVector(),
                        __indigox_pickling_version);
}

py::tuple PickleBond(const Bond bond) {
  return py::make_tuple(bond->GetAromaticity(),
                        bond->GetTag(),
//                        bond->GetMolecule(),
                        bond->GetOrder(),
                        bond->GetSourceAtom(),
                        bond->GetStereochemistry(),
                        bond->GetTargetAtom(),
                        __indigox_pickling_version);
}

py::tuple PicklePeriodicTable(const PeriodicTable table) {
  return py::make_tuple(bool(table));
}

Atom __unpickle_atom_version_0(py::tuple& t);
Bond __unpickle_bond_version_0(py::tuple& t);

Atom UnpickleAtom(py::tuple& t) {
  size_t pickled_version = t[t.size() - 1].cast<size_t>();
  switch (pickled_version) {
    case 0:
      return __unpickle_atom_version_0(t);
      
    default:
      return __unpickle_atom_version_0(t);
  }
}

Bond UnpickleBond(py::tuple& t) {
  size_t pickled_version = t[t.size() - 1].cast<size_t>();
  switch (pickled_version) {
    case 0:
      return __unpickle_bond_version_0(t);
      
    default:
      return __unpickle_bond_version_0(t);
  }
}

PeriodicTable UnpicklePeriodicTable(py::tuple& t) {
  bool real = t[0].cast<bool>();
  return real ? IXPeriodicTable::GetInstance() : PeriodicTable();
}


Atom __unpickle_atom_version_0(py::tuple& t) {
  Atom a = Atom(new IXAtom());
  a->SetAromaticity(t[0].cast<bool>());
  unsigned int Z = t[1].cast<unsigned int>();
  if (!Z) a->SetElement(IXPeriodicTable::GetInstance()->GetUndefinedElement());
  else a->SetElement(Z);
  a->SetFormalCharge(t[2].cast<int>());
  a->SetImplicitCount(t[3].cast<unsigned int>());
  a->SetIndex(t[4].cast<unsigned int>());
  a->SetName(t[5].cast<std::string>());
  return a;
}

Bond __unpickle_bond_version_0(py::tuple& t) {
  Bond bond = Bond(new IXBond());
  bond->SetAromaticity(t[0].cast<bool>());
  bond->SetTag(t[1].cast<uid_t>());
  bond->SetOrder(t[2].cast<IXBond::Order>());
  bond->SetSourceAtom(t[3].cast<Atom>());
  bond->SetStereochemistry(t[4].cast<IXBond::Stereo>());
  bond->SetTargetAtom(t[5].cast<Atom>());
  return bond;
}
