/*! \file counter.hpp */

#ifndef INDIGOX_UTILS_COUNTER_HPP
#define INDIGOX_UTILS_COUNTER_HPP

#include "numerics.hpp"

namespace indigox::utils {
  
  /*! \class IXCountableObject indigox/utils/counter.hpp
   *  \brief Template class for counting instances.
   *  \tparam T class requiring counting.
   *  \details Provides a means to count number of instances of a type created.
   *  Classes should derive from this template and call the constructor in class
   *  initalisation list. */
  template <class T>
  class IXCountableObject {
    
  public:
    /*! \brief Default constructor.
     *  \details Sets the unique id associated with the class instance.
     */
    IXCountableObject() : _id(++_count) { }
    
    /*! \brief Obtain the class instance unique id.
     *  \returns the unique id. */
    inline uid_ GetUniqueID() const { return _id; }
    
    /*! \brief Get the current instance count.
     *  \return the current number of created instances of the type. */
    static uid_ GetCurrentCount() { return _count; }
    
  private:
    //! \brief Unique ID.
    const uid_ _id;
    //! \brief Count.
    static uid_ _count;
  };
  
  template <class T>
  uid_ IXCountableObject<T>::_count = 0;
  
}

#endif /* INDIGOX_UTILS_COUNTER_HPP */
