#include <indigox/utils/modifable_object.hpp>
#include <indigox/utils/serialise.hpp>

namespace indigox::utils {
  template <typename Archive>
  void ModifiableObject::save(Archive &archive, const uint32_t) const
  {
    archive(INDIGOX_SERIAL_NVP("state", _state),
            INDIGOX_SERIAL_NVP("frozen", _frozen));
  }
  
  template <typename Archive>
  void ModifiableObject::load(Archive &archive, const uint32_t)  {
    archive(INDIGOX_SERIAL_NVP("state", _state),
            INDIGOX_SERIAL_NVP("frozen", _frozen));
  }
  INDIGOX_SERIALISE_SPLIT(ModifiableObject);
}
