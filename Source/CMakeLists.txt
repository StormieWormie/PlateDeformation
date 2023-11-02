cmake_minimum_required(VERSION 3.16.3)
project(Plate_model)

find_package(pybind11 CONFIG)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(experiment main.cpp)
add_executable(fdm fdm.cpp)
add_executable(fem fem.cpp)


pybind11_add_module(FiniteDifferenceMethod python_libs.cpp)
pybind11_add_module(FiniteElementMethod python_libs.cpp)
pybind11_add_module(Convergence python_libs.cpp)


add_custom_target(libs DEPENDS FiniteDifferenceMethod FiniteElementMethod)
add_custom_target(exe DEPENDS experiment)