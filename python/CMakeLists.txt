cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()


#install pybind11 via (from conda)
#>>pip install pybind11 
#More info here https://pybind11.readthedocs.io/en/stable/installing.html
find_package(pybind11 REQUIRED)

add_library(PyCEPS MODULE Parameterize.cpp bindings.cpp)
if(WIN32)
    set_target_properties(PyCEPS PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
    set_target_properties(PyCEPS PROPERTIES PREFIX "" SUFFIX ".so")
endif()

target_include_directories(PyCEPS PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
# borrow args.hxx directly from polyscope
target_include_directories(PyCEPS PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/polyscope/deps/args/")
target_include_directories(PyCEPS PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../src/")
target_link_libraries(PyCEPS polyscope geometry-central CEPS pybind11::module)

set_target_properties(PyCEPS PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}"
    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}"
)

#option(SPHERICAL_UNIFORMIZATION "Build spherical uniformizer" OFF)

#if(SPHERICAL_UNIFORMIZATION)
  # PETSC_DIR=/home/mark/Documents/Code/petsc PETSC_ARCH=arch-linux2-c-opt all
#  add_subdirectory(SphericalUniformization)
#endif()
