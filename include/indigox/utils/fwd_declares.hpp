// This file contains forward declarations of all classes used in indigox as
// well as their shared_ptr and weak_ptr counterparts, as needed.
#include <memory>


// Serialisation related stuff, using the cereal library
namespace cereal {
  class access;
  class PortableBinaryInputArchive;
  class PortableBinaryOutputArchive;
  class JSONInputArchive;
  class JSONOutputArchive;
  template <class T>
  class construct;
}

namespace indigox {
  
  // Molecule related
  class Molecule;
  using sMolecule = std::shared_ptr<Molecule>;
  using wMolecule = std::weak_ptr<Molecule>;
  class Atom;
  using sAtom = std::shared_ptr<Atom>;
  using wAtom = std::weak_ptr<Atom>;
  class Bond;
  using sBond = std::shared_ptr<Bond>;
  using wBond = std::weak_ptr<Bond>;
  class Angle;
  using sAngle = std::shared_ptr<Angle>;
  using wAngle = std::weak_ptr<Angle>;
  class Dihedral;
  using sDihedral = std::shared_ptr<Dihedral>;
  using wDihedral = std::weak_ptr<Dihedral>;
  class Element;
  class PeriodicTable;
  
  // CherryPicker Related
  class IXParamMolecule;
  using ParamMolecule = std::shared_ptr<IXParamMolecule>;
  using _ParamMolecule = std::weak_ptr<IXParamMolecule>;
  
  class IXParamAtom;
  using ParamAtom = std::shared_ptr<IXParamAtom>;
  using _ParamAtom = std::weak_ptr<IXParamAtom>;
  
  class IXParamBond;
  using ParamBond = std::shared_ptr<IXParamBond>;
  using _ParamBond = std::weak_ptr<IXParamBond>;
  
  class IXParamAngle;
  using ParamAngle = std::shared_ptr<IXParamAngle>;
  using _ParamAngle = std::weak_ptr<IXParamAngle>;
  
  class IXParamDihedral;
  using ParamDihedral = std::shared_ptr<IXParamDihedral>;
  using _ParamDihedral = std::weak_ptr<IXParamDihedral>;
  
  // Forcefield related
  class Forcefield;
  using sForcefield = std::shared_ptr<Forcefield>;
  using wForcefield = std::weak_ptr<Forcefield>;
  class FFAtom;
  using sFFAtom = std::shared_ptr<FFAtom>;
  using wFFAtom = std::weak_ptr<FFAtom>;
  class FFBond;
  using sFFBond = std::shared_ptr<FFBond>;
  using wFFBond = std::weak_ptr<FFBond>;
  class FFAngle;
  using sFFAngle = std::shared_ptr<FFAngle>;
  using wFFAngle = std::weak_ptr<FFAngle>;
  class FFDihedral;
  using sFFDihedral = std::shared_ptr<FFDihedral>;
  using wFFDihedral = std::weak_ptr<FFDihedral>;
  
  // Athenaeum related
  class Fragment;
  using sFragment = std::shared_ptr<Fragment>;
  using wFragment = std::weak_ptr<Fragment>;
  
  class Athenaeum;
  using sAthenaeum = std::shared_ptr<Athenaeum>;
  using wAthenaeum = std::weak_ptr<Athenaeum>;
  
  namespace algorithm {
    template <class VertType>
    struct Path;
    template <class EdgeType>
    struct EdgePath;
    
    template <class VertType>
    struct Cycle;
    template <class EdgeType>
    struct EdgeCycle;
    
    class IXCherryPicker;
    using CherryPicker = std::shared_ptr<IXCherryPicker>;
    using _CherryPicker = std::weak_ptr<IXCherryPicker>;
    
    class IXElectronAssigner;
    using ElectronAssigner = std::shared_ptr<IXElectronAssigner>;
    using _ElectronAssigner = std::weak_ptr<IXElectronAssigner>;
  }
  
  namespace graph {
    // AssignmentGraph
    class IXAssignmentGraph;
    using AssignmentGraph = std::shared_ptr<IXAssignmentGraph>;
    using _AssignmentGraph = std::weak_ptr<IXAssignmentGraph>;
    
    class IXAGVertex;
    using AGVertex = std::shared_ptr<IXAGVertex>;
    using _AGVertex = std::weak_ptr<IXAGVertex>;
    
    // MolecularGraph
    class MolecularGraph;
    using sMolecularGraph = std::shared_ptr<MolecularGraph>;
    using wMolecularGraph = std::weak_ptr<MolecularGraph>;
    class MGVertex;
    class MGEdge;
    
    // CondensedMolecularGraph
    class CondensedMolecularGraph;
    using sCondensedMolecularGraph = std::shared_ptr<CondensedMolecularGraph>;
    using wCondensedMolecularGraph = std::weak_ptr<CondensedMolecularGraph>;
    class CMGVertex;
    class CMGEdge;
  }
  
  namespace test {
    struct TestForcefield;
    struct TestFFAtom;
    struct TestFFBond;
    struct TestFFAngle;
    struct TestFFDihedral;
    
    struct TestMolecule;
    struct TestAtom;
    struct TestBond;
    struct TestAngle;
    struct TestDihedral;
    struct TestElement;
    struct TestPeriodicTable;
    
    struct TestAssignmentGraph;
    struct TestCondensedMolecularGraph;
    struct TestCondensedVertex;
    struct TestCondensedEdge;
    struct TestMolecularGraph;
    struct TestMolecularVertex;
    struct TestMolecularEdge;
  }
}
