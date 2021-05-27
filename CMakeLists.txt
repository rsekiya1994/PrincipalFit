cmake_minimum_required(VERSION 3.0)
project(PrincipalFit VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)

set(EIGEN_DIR)
find_package(Eigen3 REQUIRED)
# message(STATUS "${EIGEN3_INCLUDE_DIR}")

find_package(ROOT REQUIRED)
if(ROOT_FOUND)
    message(STATUS "ROOT FOUND : ROOT include directory = ${ROOT_INCLUDE_DIRS}")
elseif()
    message(STATUS "ROOT NOT FOUND")
endif()

if(DEFINED EIGEN_DIR)
    message(STATUS "Specified Eigen Directory=${EIGEN_DIR}")
elseif(EIGEN3_INCLUDE_DIR)
    set(EIGEN_DIR ${EIGEN3_INCLUDE_DIR})
    message(STATUS "Found Eigen : Directory=${EIGEN_DIR}")
else()
    message(STATUS "Unknown Eigen3Directories")
endif()

if(APPLE)
    set(libprefix "dylib")
else()
    set(libprefix "so")
endif()

# include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${EIGEN_DIR})
# file(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cxx)
# file(GLOB headers ${CMAKE_SOURCE_DIR}/include/*.hh)

set(libname "pcafit")
add_library(${libname} SHARED ${CMAKE_SOURCE_DIR}/src/PrincipalFit.cxx)
add_library(pcafitmacro SHARED ${CMAKE_SOURCE_DIR}/src/pca_fit_impl.cxx)
# target_include_directories(${libname} PUBLIC ${EIGEN3_INCLUDE_DIR})
target_include_directories(${libname} PUBLIC ${EIGEN_DIR})
target_include_directories(${libname} PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(${libname} PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(${libname} ${ROOT_LIBRARIES})

target_include_directories(pcafitmacro PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(pcafitmacro PUBLIC ${CMAKE_SOURCE_DIR}/macro)
target_include_directories(pcafitmacro PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(pcafitmacro ${ROOT_LIBRARIES} ${libname})

# file(CREATE_LINK "${PROJECT_BINARY_DIR}/lib${libname}.${libprefix}" ${PROJECT_SOURCE_DIR}/macro/pcafit SYMBOLIC)
# file(CREATE_LINK "${PROJECT_BINARY_DIR}/libpcafitmacro.${libprefix}" ${PROJECT_SOURCE_DIR}/macro/pcafitmac SYMBOLIC)
# add_custom_command(link_target ALL COMMAND ${CMAKE_COMMAND} -E ${PROJECT_SOURCE_DIR}/build ${PROJECT_SOURCE_DIR}/build)
# file(WRITE ${PROJECT_SOURCE_DIR}/macro/libpath.txt "lib ${PROJECT_BINARY_DIR}/lib${libname}.${libprefix}\n")
# file(APPEND ${PROJECT_SOURCE_DIR}/macro/libpath.txt "lib ${PROJECT_BINARY_DIR}/libpcafitmacro.${libprefix}\n")
# file(APPEND ${PROJECT_SOURCE_DIR}/macro/libpath.txt "inc ${PROJECT_SOURCE_DIR}/Principal.hh")

file(WRITE ${PROJECT_SOURCE_DIR}/macro/load.C "#pragma once\n")
file(APPEND ${PROJECT_SOURCE_DIR}/macro/load.C "#include <TROOT.h>\n\n")
file(APPEND ${PROJECT_SOURCE_DIR}/macro/load.C "void load() {\n")
file(APPEND ${PROJECT_SOURCE_DIR}/macro/load.C "    gROOT->ProcessLine(\"gSystem->Load(\\\"${PROJECT_BINARY_DIR}/lib${libname}.${libprefix}\\\")\");\n")
file(APPEND ${PROJECT_SOURCE_DIR}/macro/load.C "    gROOT->ProcessLine(\"gSystem->Load(\\\"${PROJECT_BINARY_DIR}/libpcafitmacro.${libprefix}\\\")\");\n")
file(APPEND ${PROJECT_SOURCE_DIR}/macro/load.C "    gROOT->ProcessLine(\"#include \\\"${PROJECT_SOURCE_DIR}/macro/pca_fit.C\\\" \");\n")
file(APPEND ${PROJECT_SOURCE_DIR}/macro/load.C "}")




