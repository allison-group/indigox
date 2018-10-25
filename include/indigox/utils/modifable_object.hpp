#ifndef modifable_object_hpp
#define modifable_object_hpp

#include <cstdint>
#include <stdexcept>

#include "fwd_declares.hpp"

namespace indigox::utils {
  
  class ModifiableObject {
    friend class cereal::access;
  public:
    using State = uint64_t;
    ModifiableObject() : _state(0), _frozen(false) { }
    inline State GetCurrentState() const { return _state; }
    inline void ModificationMade() {
      if (_frozen)
        throw std::runtime_error("Attempting to modify frozen object");
      ++_state;
    }
    inline void FreezeModifications() { _frozen = true; }
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t) const;
    template <typename Archive>
    void load(Archive& archive, const uint32_t);
    
  private:
    State _state;
    bool _frozen;
  };
  
}

#endif /* modifable_object_hpp */
