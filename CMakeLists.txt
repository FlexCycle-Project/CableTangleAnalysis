cmake_minimum_required(VERSION 3.16)
project(CableSimMin VERSION 0.1.0 LANGUAGES C CXX)

set(PROJECT_PATH ${CMAKE_SOURCE_DIR})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(EIGEN3_INCLUDE_DIR
    ${CMAKE_SOURCE_DIR}/extern/eigen
    CACHE PATH "Bundled Eigen include directory")

include_directories(
    ${PROJECT_PATH}
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/extern/GenericParameters/include
)

add_subdirectory(PositionBasedDynamics)
add_subdirectory(Utils)
add_subdirectory(Simulation)
add_subdirectory(cablegen)
