cmake_minimum_required(VERSION 3.5)
project(Calculate_CU)

add_compile_options("-g")

find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
find_package(pugixml REQUIRED)
include_directories(${MPI_INCLUDE_PATH} ${pugixml_INCLUDE_DIRS})

pkg_check_modules(PETSC_PKG REQUIRED IMPORTED_TARGET PETSc)
#pkg_check_modules(SLEPC_PKG REQUIRED IMPORTED_TARGET SLEPc)

set(CMAKE_THREAD_PREFER_PTHREAD ON)                                                                                                                                                                                                           
set(THREADS_PREFER_PTHREAD_FLAG ON)                                                                                                                                                                                                           
find_package(Threads REQUIRED) 


set(SRC_DIR src/)
set(MAIN_CPP ${SRC_DIR}main.cpp)
set(SOURCES
    ${SRC_DIR}POD_ODE_Solver.hpp
    ${SRC_DIR}POD_ODE_Solver.cpp
    ${SRC_DIR}Matrix.cpp
    ${SRC_DIR}Matrix.hpp
)
set(XML_FILE ../POD_Para.xml)
set(ALL_SOURCES
    ${MAIN_CPP}
    ${SOURCES}
)

add_compile_options("-lpthread -fpermissive")
set(GCC_CXX_FLAGS ${GCC_CXX_FLAGS} "-std=c++11 -lpthread -fpermissive")
### MAIN PROGRAMS #########################################

add_executable(${PROJECT_NAME} ${ALL_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}>
    ${MPI_INCLUDE_PATH}
    ${pugixml_INCLUDE_DIRS}
)
target_link_libraries(${PROJECT_NAME} PRIVATE 
    Threads::Threads
    ${MPI_LIBRARIES}
    PkgConfig::PETSC_PKG
    ${pugixml_LIBRARIES}
    
    #PkgConfig::SLEPC_PKG
)
target_link_libraries(${PROJECT_NAME} PRIVATE /usr/lib/x86_64-linux-gnu/libpugixml.a)



