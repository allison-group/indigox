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
#OPTION(BUILD_EXAMPLES "Build C++ example targets" ON) # Currently not used - examples always built
OPTION(INSTALL_BOOST "Install boost headers if using provided" OFF)

# Default to release mode build. Does Debug actually do anything?
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
  MESSAGE(STATUS "Defaulting to Release mode build. Use -DCMAKE_BUILD_TYPE=Debug to build in debug mode, though it doesn't really do anything currently.")
ENDIF()

# Setup some compile flags, for all versions and specific versions
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -mavx2")
IF("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF()
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -O0")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/)

# Set default install prefix to /usr/local/opt/indigox
SET(CMAKE_INSTALL_PREFIX /usr/local/opt/${PROJECT_NAME})

# Setup install directories
INCLUDE(GNUInstallDirs)
ADD_DEFINITIONS(-DIX_DATA_DIRECTORY="${CMAKE_INSTALL_FULL_DATAROOTDIR}")
ADD_DEFINITIONS(-DINDIGOX_DISABLE_SANITY_CHECKS=1)
ADD_DEFINITIONS(-DINDIGOX_DISABLE_TESTS=1)

# List all the files to compile into the indigox library
SET(INDIGOX_LIB_SRCS
    src/algorithm/cherrypicker.cpp
    src/algorithm/graph/connectivity.cpp
    src/algorithm/graph/cycles.cpp
    src/algorithm/graph/isomorphism.cpp
    src/algorithm/graph/paths.cpp
    src/classes/angle.cpp
    src/classes/athenaeum.cpp
    src/classes/atom.cpp
    src/classes/bond.cpp
    src/classes/dihedral.cpp
    src/classes/forcefield.cpp
    src/classes/molecule.cpp
    src/classes/parameterised.cpp
    src/classes/periodictable.cpp
    src/classes/residue.cpp
    src/graph/condensed.cpp
    src/graph/molecular.cpp
    src/utils/common.cpp
    )

# Build the library
ADD_LIBRARY(indigox STATIC ${INDIGOX_LIB_SRCS})

TARGET_COMPILE_FEATURES(indigox PRIVATE cxx_std_17)
TARGET_COMPILE_FEATURES(indigox PRIVATE cxx_inline_namespaces)

# Add external paths to include directories
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/EASTL/include)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/cereal/include)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/doctest)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/eigen)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/RI)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/json/single_include)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/indigo-bondorder/include)

# We provide a boost version so only use system boost if it is newer than
# provided boost
# Sometime I need to figure out a minimum required version for it. :S
FIND_PACKAGE(Boost 1.66.0)
IF(NOT Boost_FOUND)
  MESSAGE(STATUS "Using provided Boost version.")
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/boost/include)
  IF (INSTALL_BOOST)
    INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/external/boost/include/boost
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  ENDIF()
ELSE()
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

# Build bindings if requested
IF(BUILD_PYTHON)
  ADD_DEFINITIONS(-DINDIGOX_DISABLE_TESTS=1)
  ADD_SUBDIRECTORY(src/python)
ELSE()
  MESSAGE(STATUS "Building of Python bindings disabled. Your loss.")
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

# Install the indigo-bondorder library too
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/external/indigo-bondorder)
TARGET_LINK_LIBRARIES(indigox indigo-bondorder)

ADD_EXECUTABLE(nitrobenzoate_example examples/CherryPicker/nitrobenzoate_example.cpp)
TARGET_LINK_LIBRARIES(nitrobenzoate_example indigox)
TARGET_LINK_LIBRARIES(nitrobenzoate_example stdc++fs) # needed for std:filesystem

ADD_EXECUTABLE(general_example examples/CherryPicker/general_example.cpp)
TARGET_LINK_LIBRARIES(general_example indigox)
TARGET_LINK_LIBRARIES(general_example stdc++fs) # needed for std:filesystem

FILE(COPY data DESTINATION .)

# Install library
INSTALL(TARGETS indigox DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install examples
INSTALL(TARGETS nitrobenzoate_example DESTINATION ${CMAKE_INSTALL_LIBDIR})
INSTALL(TARGETS general_example DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install header files. This should includes boost if the repo version is being used
#INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/indigox
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install data files
INSTALL(FILES data/periodictable.json
       DESTINATION ${CMAKE_INSTALL_DATAROOTDIR})
# Install documentation
#INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION ${CMAKE_INSTALL_DOCDIR} OPTIONAL)

# Uninstall target
IF(NOT TARGET uninstall)
  CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                  IMMEDIATE @ONLY)
  ADD_CUSTOM_TARGET(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake VERBATIM)
ENDIF()
