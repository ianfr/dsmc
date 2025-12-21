cmake_minimum_required(VERSION 3.24)
project(dsmc)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

include_directories(
    depends
)

add_executable(dsmc main.cpp Cell.cpp Cell.h Particle.h Grid.cpp Grid.h)

# https://stackoverflow.com/a/52230189
target_link_libraries (dsmc Eigen3::Eigen)
