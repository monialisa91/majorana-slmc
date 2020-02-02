cmake_minimum_required(VERSION 3.15)
project(rashba)

set(CMAKE_CXX_STANDARD 14)

add_executable(rashba main.cpp hamiltonian.h)
target_link_libraries(rashba armadillo)
