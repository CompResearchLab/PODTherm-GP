cmake_minimum_required(VERSION 3.5)
project(Therm_FEM)
#project(Xeon_CPU)
cmake_policy(SET CMP0004 NEW)

#set(CMAKE_MODULE_PATH "/home/tony/other_research/fem_fenics/fem_cpp/cmake/cmake/modules")


find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

set(SRC_DIR src/)
set(MAIN_CPP ${SRC_DIR}main.cpp)
#include_directories(./usr/include)

#set(MAIN_CPP ${SRC_DIR} mai)
set(OTHER_SOURCES
    ${SRC_DIR}helpers.hpp
    ${SRC_DIR}helpers.cpp
    ${SRC_DIR}Space.h
    #${SRC_DIR}expressions_AMD_finesource.hpp
    ${SRC_DIR}expressions.hpp
    ${SRC_DIR}boundary.hpp
)

set(ALL_SOURCES
    ${MAIN_CPP}
    ${OTHER_SOURCES}
)

add_executable(${PROJECT_NAME} ${ALL_SOURCES})

#set(CMAKE_CPP_OPTIONS -mfpmath=sse)
#add_compiler_options(-mno-fancy-math-387)
#set(CMAKE_CPP_FLAGS -mno-fancy-math-387)
target_link_libraries(${PROJECT_NAME} dolfin)
target_link_libraries(${PROJECT_NAME} /usr/lib/x86_64-linux-gnu/libpugixml.a)

target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<BUILD_INTERFACE:/usr/lib/petsc/include>
    $<BUILD_INTERFACE:/usr/lib/slepc/include>
)
