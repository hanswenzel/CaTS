cmake_minimum_required(VERSION 3.16...3.21)

if(WITH_G4OPTICKS AND WITH_CXG4OPTICKS)
  message(FATAL_ERROR "CaTS cmake illegal set of parameters. One can't use opticks and legacy opticks at the same time. Either -DWITH_G4OPTICKS or -DWITH_CXG4OPTICKS must be set to off and can't be set on at the same time.
the following combinations are allowed:
-DWITH_G4OPTICKS=OFF -DWITH_CXG4OPTICKS=OFF : no opticks just Geant4 optical physics
-DWITH_G4OPTICKS=ON  -DWITH_CXG4OPTICKS=OFF : build using legacy opticks based on NVIDIA OptiX 6
-DWITH_G4OPTICKS=OFF -DWITH_CXG4OPTICKS=ON  : build using current opticks based on NVIDIA OptiX > 7
")
endif()

set(name CaTS)
project(${name} VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
option(WITH_ROOT "Build example with ROOT" ON)
option(WITH_G4OPTICKS "Build example with OPTICKS" OFF)
option(WITH_CXG4OPTICKS "Build example with OPTICKS" OFF)
option(G4ANALYSIS_USE "Build example with Analysis" ON)

# enable_testing()
include(CTest)

# ----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# ----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

# ----------------------------------------------------------------------------
# CaTS requires shared libraries
#
if(NOT Geant4_shared_FOUND)
  message(FATAL_ERROR "CaTS must use shared libraries")
endif()

if(WITH_G4OPTICKS)
  message(STATUS "WITH_G4OPTICKS is set")
  list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
  include(OpticksBuildOptions)
  find_package(G4OK CONFIG REQUIRED)
  include(CaTSCXXFlags)
elseif(WITH_CXG4OPTICKS)
  message(STATUS "WITH_CXG4OPTICKS is set")
  list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
  message(STATUS "*************************************************************")
  message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
  message(STATUS "OPTICKS_HOME:  $ENV{OPTICKS_HOME}")
  message(STATUS "*************************************************************")
  include(OpticksBuildOptions)
  find_package(G4CX REQUIRED CONFIG)
  find_package(U4 REQUIRED CONFIG)
  find_package(ExtG4 REQUIRED CONFIG)
  find_package(CSG_GGeo REQUIRED CONFIG)

  # include(OpticksBuildOptions)
  include(CaTSCXXFlags)

# find_package(CSGOptiX REQUIRED CONFIG)
else()
  message(STATUS "CaTS is built without opticks support")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
  include(CaTSCXXFlags)
endif()

if(WITH_ROOT)
  find_package(ROOT COMPONENTS EG REQUIRED)

  # ROOT version 6 required
  if(ROOT_FOUND)
    STRING(REGEX MATCH "6.*" VERSION6MATCH ${ROOT_VERSION})

    if(NOT VERSION6MATCH)
      message(FATAL_ERROR "${name} requires ROOT 6")
    endif()
  endif()

  # Include ROOT's CMake functions for dictionary generation
  # since root6.20, the file is renamed and included by default, so include
  # only when we find the *old* name
  if(EXISTS "${ROOT_DIR}/modules/RootNewMacros.cmake")
    include("${ROOT_DIR}/modules/RootNewMacros.cmake")
  endif()
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

if(WITH_ROOT)
  # ----------------------------------------------------------------------------
  # Locate sources and headers for this project
  #
  include_directories(${PROJECT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIRS})

  # ----------------------------------------------------------------------------
  # Generate dictionaries, add ROOT libraries properties
  #
  # message(STATUS "root libraries: " ${ROOT_LIBRARIES})
  REFLEX_GENERATE_DICTIONARY(CaTSClasses include/CaTSClasses.hh SELECTION xml/selection.xml)
  add_library(CaTSClassesDict SHARED ${sources} CaTSClasses.cxx)
  set(libsuffix .so)
  set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES} SUFFIX ${libsuffix})
  set_target_properties(CaTSClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})
  target_link_libraries(CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
else()
  include_directories(${PROJECT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIR})
endif()

file(GLOB detectors ${PROJECT_SOURCE_DIR}/gdml/*.gdml)
file(GLOB fields ${PROJECT_SOURCE_DIR}/gdml/*.root)
file(GLOB schemas ${PROJECT_SOURCE_DIR}/gdml/*.xsd)
file(GLOB scripts ${PROJECT_SOURCE_DIR}/scripts/*)
file(GLOB macros ${PROJECT_SOURCE_DIR}/macros/*.mac)
file(GLOB runscripts ${PROJECT_SOURCE_DIR}/scripts/run.sh ${PROJECT_SOURCE_DIR}/scripts/check.sh)

if(WITH_G4OPTICKS)
  add_executable(${name} ${name}.cc ${sources} ${headers})
  target_compile_definitions(${name} PRIVATE WITH_G4OPTICKS WITH_ROOT)
  target_link_libraries(${name} CaTSClassesDict Opticks::G4OK ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_test(${name}Test ${name} -g simpleLArTPC.gdml -m test_G4.mac)
elseif(WITH_CXG4OPTICKS)
  add_executable(${name} ${name}.cc ${sources} ${headers})
  target_compile_definitions(${name} PRIVATE WITH_CXG4OPTICKS WITH_ROOT)
  target_link_libraries(${name} CaTSClassesDict Opticks::G4CX Opticks::U4 Opticks::ExtG4 Opticks::CSG_GGeo ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_test(${name}Test ${name} -g simpleLArTPC.gdml -m test_G4.mac)
else()
  add_executable(${name} ${name}.cc ${sources} ${headers})
  target_compile_definitions(${name} PRIVATE WITH_ROOT)
  target_link_libraries(${name} CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_test(${name}Test ${name} -g simpleLArTPC.gdml -m test_G4.mac)
endif()

if(WITH_ROOT)
  add_executable(readPhotonHits readPhotonHits.cc ${sources} ${headers})
  target_link_libraries(readPhotonHits CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_executable(readCalorimeterHits readCalorimeterHits.cc ${sources} ${headers})
  target_link_libraries(readCalorimeterHits CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_executable(readMscHits readMscHits.cc ${sources} ${headers})
  target_link_libraries(readMscHits CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_executable(readDRCalorimeterHits readDRCalorimeterHits.cc ${sources} ${headers})
  target_link_libraries(readDRCalorimeterHits CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_executable(readTrackerHits readTrackerHits.cc ${sources} ${headers})
  target_link_libraries(readTrackerHits CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_executable(readSimTrajectory readSimTrajectory.cc ${sources} ${headers})
  target_link_libraries(readSimTrajectory CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_executable(readSimEnergyDeposit readSimEnergyDeposit.cc ${sources} ${headers})
  target_link_libraries(readSimEnergyDeposit CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  add_executable(readInteractionHits readInteractionHits.cc ${sources} ${headers})
  target_link_libraries(readInteractionHits CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
  link_directories(${ROOT_LIBRARY_DIR})
endif()

install(TARGETS ${name} DESTINATION bin)
install(FILES ${detectors} ${fields} ${schemas} ${scripts} ${macros} DESTINATION bin)
install(FILES ${detectors} DESTINATION bin/gdml)
install(FILES ${macros} DESTINATION bin/macros)
install(PROGRAMS ${runscripts} DESTINATION bin)

if(WITH_ROOT)
  install(TARGETS CaTSClassesDict DESTINATION bin)
  install(TARGETS readPhotonHits DESTINATION bin)
  install(TARGETS readCalorimeterHits DESTINATION bin)
  install(TARGETS readDRCalorimeterHits DESTINATION bin)
  install(TARGETS readSimTrajectory DESTINATION bin)
  install(TARGETS readSimEnergyDeposit DESTINATION bin)
  install(TARGETS readMscHits DESTINATION bin)
  install(TARGETS readTrackerHits DESTINATION bin)
  install(TARGETS readInteractionHits DESTINATION bin)
  install(FILES ${PROJECT_BINARY_DIR}/libCaTSClassesDict.so DESTINATION bin)
  install(FILES ${PROJECT_BINARY_DIR}/CaTSClasses_rdict.pcm DESTINATION bin)
endif()

message(STATUS "G4ANALYSIS_USE:   ${G4ANALYSIS_USE}")
message(STATUS "WITH_G4OPTICKS:   ${WITH_G4OPTICKS}")
message(STATUS "WITH_CXG4OPTICKS: ${WITH_CXG4OPTICKS}")
message(STATUS "WITH_ROOT:        ${WITH_ROOT}")
