cmake_minimum_required(VERSION 3.21)
project(Geant4Wrap)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

#---Find JlCxx package-------------------------------------------------------------
find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")
#---Find Geant4---------------------------------------------------------------------
find_package(Geant4 REQUIRED gdml)
include(${Geant4_USE_FILE})

#  Files are in gen folder
file(REAL_PATH ${CMAKE_SOURCE_DIR}/src SOURCE_DIR)
file(GLOB GEN_SOURCES CONFIGURE_DEPENDS  ${SOURCE_DIR}/cpp/Jl*.cxx)
add_library(Geant4Wrap SHARED ${SOURCE_DIR}/cpp/Geant4Wrap.cxx ${SOURCE_DIR}/cpp/jlGeant4.cxx ${GEN_SOURCES})

#target_compile_definitions(Geant4Wrap PUBLIC VERBOSE_IMPORT)
target_include_directories(Geant4Wrap PRIVATE ${SOURCE_DIR})
target_link_libraries(Geant4Wrap JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl ${Geant4_LIBRARIES})

install(TARGETS Geant4Wrap
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION lib)
install(FILES src/jl/Geant4-export.jl DESTINATION .)
