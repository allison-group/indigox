#ifndef INDIGOX_UTILS_SERIALISE_HPP
#define INDIGOX_UTILS_SERIALISE_HPP

#define CEREAL_SERIALIZE_FUNCTION_NAME Serialise

#include <type_traits>

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/bitset.hpp>
#include <cereal/types/eastl_bitset.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

#define INDIGOX_SERIALISE(Class_name) \
template void Class_name::Serialise<cereal::PortableBinaryInputArchive>\
(cereal::PortableBinaryInputArchive&, const uint32_t); \
template void Class_name::Serialise<cereal::PortableBinaryOutputArchive>\
(cereal::PortableBinaryOutputArchive&, const uint32_t); \
template void Class_name::Serialise<cereal::JSONInputArchive>\
(cereal::JSONInputArchive&, const uint32_t); \
template void Class_name::Serialise<cereal::JSONOutputArchive>\
(cereal::JSONOutputArchive&, const uint32_t);

#define INDIGOX_SERIAL_NVP(name, value) cereal::make_nvp<Archive>(name, value)

#define INDIGOX_IS_OUTPUT_ARCHIVE \
std::is_base_of<cereal::detail::OutputArchiveBase, Archive>::value

#define INDIGOX_IS_INPUT_ARCHIVE \
std::is_base_of<cereal::detail::InputArchiveBase, Archive>::value

#endif /* INDIGOX_UTILS_SERIALISE_HPP */
