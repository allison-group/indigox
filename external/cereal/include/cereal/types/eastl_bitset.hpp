/*! \file eastl_bitset.hpp
    \brief Support for types found in \<EASTL/bitset.h\>
    \ingroup STLSupport */
#ifndef CEREAL_TYPES_EASTL_BITSET_HPP_
#define CEREAL_TYPES_EASTL_BITSET_HPP_

#include "cereal/cereal.hpp"
#include <EASTL/bitset.h>

namespace cereal
{
  //! Serializing (save) for eastl::bitset
  template <class Archive, size_t N, class W> inline
  void CEREAL_SAVE_FUNCTION_NAME( Archive & ar, eastl::bitset<N, W> const & bits )
  {
    using Bitset_t = eastl::bitset<N, W>;
    W words[Bitset_t::kWordCount];
    for (size_t i = 0; i < Bitset_t::kWordCount; ++i) words[i] = bits.data()[i];
    ar( CEREAL_NVP_("data", words) );
  }
  
  //! Serializing (load) for eastl::bitset
  template <class Archive, size_t N, class W> inline
  void CEREAL_LOAD_FUNCTION_NAME( Archive & ar, eastl::bitset<N, W> & bits )
  {
    using Bitset_t = eastl::bitset<N, W>;
    W words[Bitset_t::kWordCount];
    ar( CEREAL_NVP_("data", words) );
    for (size_t i = 0; i < Bitset_t::kWordCount; ++i) bits.data()[i] = words[i];
  }
} // namespace cereal

#endif // CEREAL_TYPES_EASTL_BITSET_HPP_
