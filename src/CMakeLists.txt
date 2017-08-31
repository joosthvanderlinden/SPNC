# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( spnc )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package

# Creating entries for target: spnc
# ############################

add_executable( spnc main.cpp globals.h 
                configure.cpp configure.h
                pack.cpp pack.h 
                transform.cpp transform.h 
                intersect.cpp intersect.h 
                disjointset.cpp disjointset.h 
                particles.cpp particles.h 
                porenet.cpp porenet.h 
                contactnet.cpp contactnet.h 
                triangulate.cpp triangulate.h 
                rangetree.cpp rangetree.h)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS spnc)

# Link the executable to CGAL and third-party libraries
target_link_libraries(spnc   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

# Set the required C++ standard automatically 
# (uncomment this if you use CMake 3.1.0 or later and/or your compiler is not set to C++11 by default)
# target_compile_features(spnc PRIVATE cxx_range_for)