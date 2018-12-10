/*! \file eigen.hpp
    \brief Support for types found in \<Eigen/Dense\>
    \ingroup STLSupport */
#ifndef CEREAL_TYPES_EIGEN_HPP_
#define CEREAL_TYPES_EIGEN_HPP_

#include "cereal/cereal.hpp"
#include <Eigen/Dense>

namespace cereal
{
  //! Serializing (save) for Eigen::Vector3d
  template <class Archive> inline
  void CEREAL_SAVE_FUNCTION_NAME( Archive & ar, Eigen::Vector3d const & vec )
  {
    ar(vec[0], vec[1], vec[2]);
  }
  
  //! Serializing (load) for Eigen::Vector3d
  template <class Archive> inline
  void CEREAL_LOAD_FUNCTION_NAME( Archive & ar, Eigen::Vector3d & vec )
  {
    ar(vec[0], vec[1], vec[2]);
  }
} // namespace cereal

#endif // CEREAL_TYPES_EIGEN_HPP_
