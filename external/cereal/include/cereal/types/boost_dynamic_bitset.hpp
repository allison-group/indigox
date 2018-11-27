/*! \file boost_dynamic_bitset.hpp
    \brief Support for types found in \<boost/dynamic_bitset/dynamic_bitset.hpp\>
    \ingroup STLSupport */
#ifndef CEREAL_TYPES_BOOST_DYNAMIC_BITSET_HPP_
#define CEREAL_TYPES_BOOST_DYNAMIC_BITSET_HPP_

#include "cereal/cereal.hpp"
#include <iterator>
#include <vector>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

namespace cereal
{
  //! Serializing (save) for boost::dynamic_bitset
  template <class Archive> inline
  void CEREAL_SAVE_FUNCTION_NAME( Archive & ar, boost::dynamic_bitset<> const & bits )
  {
    using Bitset_t = boost::dynamic_bitset<>;
    using block = Bitset_t::block_type;
    std::vector<block> blocks; blocks.reserve(bits.num_blocks());
    boost::to_block_range(bits, std::back_inserter(blocks));
    ar(CEREAL_NVP_("blocks", blocks));
  }
  
  //! Serializing (load) for boost::dynamic_bitset
  template <class Archive> inline
  void CEREAL_LOAD_FUNCTION_NAME( Archive & ar, boost::dynamic_bitset<> & bits )
  {
    using Bitset_t = boost::dynamic_bitset<>;
    using block = Bitset_t::block_type;
    std::vector<block> blocks;
    ar(CEREAL_NVP_("blocks", blocks));
    bits.reserve(blocks.size() * Bitset_t::bits_per_block);
    boost::from_block_range(blocks.begin(), blocks.end(), bits);
  }
} // namespace cereal

#endif // CEREAL_TYPES_BOOST_DYNAMIC_BITSET_HPP_
