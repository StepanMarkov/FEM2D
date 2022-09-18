
cmake_minimum_required (VERSION 3.8)

#find_package(Eigen3 REQUIRED)
#set(Eigen3_DIR "${Eigen3_DIR}/include/eigen3/Eigen")
#include_directories(${Eigen3_INCLUDE_DIRS})

add_executable (FEM 
"FEM.cpp" 
"FEM.h"
"FEM_BASE.cpp"
"FEM_MESH.cpp"
"FEM_EQUATION.cpp"
"FEM_SOLVER.cpp"
"FEM_RESULTS.cpp"
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET FEM PROPERTY CXX_STANDARD 20)
endif()

# TODO: Добавьте тесты и целевые объекты, если это необходимо.
