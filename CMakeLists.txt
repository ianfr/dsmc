cmake_minimum_required(VERSION 3.24)
project(dsmc)

set(CMAKE_CXX_STANDARD 17)

add_executable(dsmc main.cpp Cell.cpp Cell.h Particle.h Grid.cpp Grid.h)
