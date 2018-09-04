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
  
  class IXMolecule;
  using Molecule = std::shared_ptr<IXMolecule>;
  using _Molecule = std::weak_ptr<IXMolecule>;
  
  namespace algorithm {
    template <class VertType>
    struct Path;
    template <class EdgeType>
    struct EdgePath;
    
    template <class VertType>
    struct Cycle;
    template <class EdgeType>
    struct EdgeCycle;
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
    
  }
}
