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
  class IXMolecule;
  using Molecule = std::shared_ptr<IXMolecule>;
  using _Molecule = std::weak_ptr<IXMolecule>;
  
  class IXAtom;
  using Atom = std::shared_ptr<IXAtom>;
  using _Atom = std::weak_ptr<IXAtom>;
  
  class IXBond;
  using Bond = std::shared_ptr<IXBond>;
  using _Bond = std::weak_ptr<IXBond>;
  
  class IXAngle;
  using Angle = std::shared_ptr<IXAngle>;
  using _Angle = std::weak_ptr<IXAngle>;
  
  class IXDihedral;
  using Dihedral = std::shared_ptr<IXDihedral>;
  using _Dihedral = std::weak_ptr<IXDihedral>;
  
  class IXElement;
  using Element = std::shared_ptr<IXElement>;
  using _Element = std::weak_ptr<IXElement>;
  
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
  class IXForcefield;
  using Forcefield = std::shared_ptr<IXForcefield>;
  using _Forcefield = std::weak_ptr<IXForcefield>;
  
  class IXFFAtom;
  using FFAtom = std::shared_ptr<IXFFAtom>;
  using _FFAtom = std::weak_ptr<IXFFAtom>;
  
  class IXFFBond;
  using FFBond = std::shared_ptr<IXFFBond>;
  using _FFBond = std::weak_ptr<IXFFBond>;
  
  class IXFFAngle;
  using FFAngle = std::shared_ptr<IXFFAngle>;
  using _FFAngle = std::weak_ptr<IXFFAngle>;
  
  class IXFFDihedral;
  using FFDihedral = std::shared_ptr<IXFFDihedral>;
  using _FFDihedral = std::weak_ptr<IXFFDihedral>;
  
  // Athenaeum related
  class IXFragment;
  using Fragment = std::shared_ptr<IXFragment>;
  using _Fragment = std::weak_ptr<IXFragment>;
  
  class IXAthenaeum;
  using Athenaeum = std::shared_ptr<IXAthenaeum>;
  using _Athenaeum = std::weak_ptr<IXAthenaeum>;
  
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
    class IXMolecularGraph;
    using MolecularGraph = std::shared_ptr<IXMolecularGraph>;
    using _MolecularGraph = std::weak_ptr<IXMolecularGraph>;
    
    class IXMGVertex;
    using MGVertex = std::shared_ptr<IXMGVertex>;
    using _MGVertex = std::weak_ptr<IXMGVertex>;
    
    class IXMGEdge;
    using MGEdge = std::shared_ptr<IXMGEdge>;
    using _MGEdge = std::weak_ptr<IXMGEdge>;
    
    // CondensedMolecularGraph
    class IXCondensedMolecularGraph;
    using CondensedMolecularGraph = std::shared_ptr<IXCondensedMolecularGraph>;
    using _CondensedMolecularGraph = std::weak_ptr<IXCondensedMolecularGraph>;
    
    class IXCMGVertex;
    using CMGVertex = std::shared_ptr<IXCMGVertex>;
    using _CMGVertex = std::weak_ptr<IXCMGVertex>;
    
    class IXCMGEdge;
    using CMGEdge = std::shared_ptr<IXCMGEdge>;
    using _CMGEdge = std::weak_ptr<IXCMGEdge>;
  }
  
  namespace test {
    struct TestForcefield;
    struct TestFFAtom;
    struct TestFFBond;
    struct TestFFAngle;
    struct TestFFDihedral;
    struct TestMolecule;
    struct TestAtom;
    struct TestAssignmentGraph;
    struct TestCondensedMolecularGraph;
    struct TestCondensedVertex;
    struct TestCondensedEdge;
  }
}
