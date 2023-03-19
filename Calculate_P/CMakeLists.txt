cmake_minimum_required(VERSION 3.5)
project(calculate_P_Matrix)
cmake_policy(SET CMP0004 NEW)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

set(SRC_DIR src/)
set(MAIN_CPP ${SRC_DIR}main.cpp)
set(OTHER_SOURCES
    ${SRC_DIR}helpers.hpp
    ${SRC_DIR}helpers.cpp
    ${SRC_DIR}Space.h
    ${SRC_DIR}expressions.hpp
    ${SRC_DIR}boundary.hpp
    ${SRC_DIR}Matrix.hpp
    ${SRC_DIR}Matrix.cpp
)

set(ALL_SOURCES
    ${MAIN_CPP}
    ${OTHER_SOURCES}
)

add_executable(${PROJECT_NAME} ${ALL_SOURCES})

set(GCC_CXX_FLAGS ${GCC_CXX_FLAGS} "-static")
target_link_libraries(${PROJECT_NAME} dolfin)
target_link_libraries(${PROJECT_NAME} /usr/lib/x86_64-linux-gnu/libpugixml.a)
target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<BUILD_INTERFACE:/usr/lib/petsc/include>
    $<BUILD_INTERFACE:/usr/lib/slepc/include>
)
