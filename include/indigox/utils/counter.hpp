/*! \file counter.hpp */

#ifndef INDIGOX_UTILS_COUNTER_HPP
#define INDIGOX_UTILS_COUNTER_HPP

#include <cstdint>

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
    inline uint64_t GetUniqueID() const { return _id; }
    
  private:
    //! \brief Unique ID.
    const uint64_t _id;
    //! \brief Count.
    static uint32_t _count;
  };
  
  template <class T>
  uint32_t IXCountableObject<T>::_count = 0;
  
}

#endif /* INDIGOX_UTILS_COUNTER_HPP */
