cmake_minimum_required( VERSION 3.25 )
project(global_optimization_algorithms)

set(CMAKE_CXX_STANDARD 20)

add_executable(global_optimization_algorithms main.cpp Solve.cpp)
