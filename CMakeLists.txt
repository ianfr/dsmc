cmake_minimum_required(VERSION 3.24)
project(dsmc)

set(CMAKE_CXX_STANDARD 17)

if (APPLE)
    include_directories(
        /opt/homebrew/Cellar/eigen/5.0.1/include
    )
endif()

add_executable(dsmc main.cpp Cell.cpp Cell.h Particle.h Grid.cpp Grid.h)
