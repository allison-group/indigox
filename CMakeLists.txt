# Should figure out minimum needed at some point. 
CMAKE_MINIMUM_REQUIRED(VERSION 3.7)
CMAKE_POLICY(VERSION 3.7)

PROJECT(indigox VERSION 0.7.2.0)

# We require c++-17 compiler
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Set defaults for available options
OPTION(BUILD_PYTHON "Build Python bindings" ON)
OPTION(BUILD_EXAMPLES "Build a few example programs" OFF)

# Default to release mode build. Does Debug actually do anything?
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
  MESSAGE(STATUS "Defaulting to Release mode build. Use -DCMAKE_BUILD_TYPE=Debug to build in debug mode, though it doesn't really do anything currently.")
ENDIF()

# Setup some compile flags, for all versions and specific versions
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-deprecated-declarations -Wno-deprecated")
IF("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF()
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -DDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/)

# Set default install prefix to /usr/local/opt/indigox
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)    # Requires Cmake >= 3.7
  # Becuase I don't like populating /usr/local, make it in /usr/local/opt by default
  SET(CMAKE_INSTALL_PREFIX "/usr/local/opt/${PROJECT_NAME}" CACHE PATH "..." FORCE)
ELSE()
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME} CACHE PATH "..." FORCE)
ENDIF()
# Setup install directories
INCLUDE(GNUInstallDirs)
ADD_DEFINITIONS(-DDFT_DATA_DIR=${CMAKE_INSTALL_FULL_DATAROOTDIR})
ADD_DEFINITIONS(-DFCSCORE=uint32_t)

# List all the files to compile into the indigox library
SET(INDIGOX_LIB_SRCS
    src/classes/angle.cpp
    src/classes/atom.cpp
    src/classes/bond.cpp
    src/classes/dihedral.cpp
    src/classes/molecule.cpp
    src/classes/periodictable.cpp
    src/graph/assignment.cpp
    src/graph/molecular.cpp
    src/utils/common.cpp
    src/utils/filereader.cpp
    src/utils/options.cpp)

# Build the library
ADD_LIBRARY(indigox STATIC ${INDIGOX_LIB_SRCS})

# We provide a boost version so only use system boost if it is newer than
# provided boost
# Sometime I need to figure out a minimum required version for it. :S
FIND_PACKAGE(Boost 1.66.0)
IF(NOT Boost_FOUND)
  MESSAGE(STATUS "Using provided Boost version.")
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/boost/include)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/external/boost/include/boost 
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
ELSE()
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

# Build bindings if requested
IF(BUILD_PYTHON)
  ADD_SUBDIRECTORY(src/python)
ELSE()
  MESSAGE(STATUS "Building of Python bindings disabled.")
ENDIF()

IF(BUILD_EXAMPLES)
    ADD_EXECUTABLE(cxx_example examples/cxx_example.cpp)
    TARGET_INCLUDE_DIRECTORIES(cxx_example BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/include)
    TARGET_LINK_LIBRARIES(cxx_example indigox)
    ADD_CUSTOM_TARGET(examples ALL)
    ADD_DEPENDENCIES(examples cxx_example)
    INSTALL(TARGETS cxx_example DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
    INSTALL(FILES examples/cxx_example.cpp examples/python_example.py
            DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
ENDIF()

# Build documentation
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  SET(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
  SET(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  CONFIGURE_FILE(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  ADD_CUSTOM_TARGET(doc 
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with doxygen"
  VERBATIM )
  MESSAGE(STATUS "Doxygen located. Use 'make doc' to build documentation.")
ELSE(DOXYGEN_FOUND)
  MESSAGE(STATUS "Unable to locate doxygen. Documentation will NOT be built.")
ENDIF(DOXYGEN_FOUND)

FILE(COPY data DESTINATION .)
FILE(COPY examples/python_example.py DESTINATION examples)

# Install examples
# Install library
INSTALL(TARGETS indigox DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install header files. This should includes boost if the repo version is being used
#INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/indigox
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install data files
#INSTALL(FILES data/atomenergies_def2svpd.txt
#        data/atomenergies_def2tzvppd.txt
#        data/bondenergies_def2svpd.txt
#        data/bondenergies_def2svpd_withbsse.txt
#        data/bondenergies_def2tzvppd.txt
#        data/bondenergies_def2tzvppd_withbsse.txt
#        data/periodictable.txt
#        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR})
# Install documentation
#INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION ${CMAKE_INSTALL_DOCDIR} OPTIONAL)

# Uninstall target
IF(NOT TARGET uninstall)
  CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                  IMMEDIATE @ONLY)
  ADD_CUSTOM_TARGET(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake VERBATIM)
ENDIF()
