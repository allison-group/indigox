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

#define __hr_out  cereal::JSONOutputArchive
#define __hr_in   cereal::JSONInputArchive
#define __cr_out  cereal::PortableBinaryOutputArchive
#define __cr_in   cereal::PortableBinaryInputArchive

#define __ix_serialise(class_name, archive_t) \
template void class_name::Serialise<archive_t> (archive_t&, const uint32_t);

#define __ix_serialise_split(class_name, in_archive_t, out_archive_t) \
template void class_name::save<out_archive_t> (out_archive_t&, const uint32_t) const; \
template void class_name::load_and_construct<in_archive_t> \
(in_archive_t&, cereal::construct<class_name>&, const uint32_t); 

#define INDIGOX_SERIALISE(class_name) \
__ix_serialise(class_name, __hr_out) __ix_serialise(class_name, __cr_out) \
__ix_serialise(class_name, __hr_in) __ix_serialise(class_name, __cr_in)

#define INDIGOX_SERIALISE_SPLIT(class_name) \
__ix_serialise_split(class_name, __hr_in, __hr_out) \
__ix_serialise_split(class_name, __cr_in, __cr_out)

#define INDIGOX_SERIAL_NVP(name, value) cereal::make_nvp<Archive>(name, value)

#define INDIGOX_SERIALISE_VERSION(class_name, version) \
CEREAL_CLASS_VERSION(class_name, version);

#define INDIGOX_IS_OUTPUT_ARCHIVE(archive_t) \
std::is_base_of<cereal::detail::OutputArchiveBase, archive_t>::value

#define INDIGOX_IS_INPUT_ARCHIVE(archive_t) \
std::is_base_of<cereal::detail::InputArchiveBase, archive_t>::value

#define INDIGOX_IS_HUMAN_READABLE(archive_t) \
cereal::traits::is_text_archive<archive_t>::value

#endif /* INDIGOX_UTILS_SERIALISE_HPP */
