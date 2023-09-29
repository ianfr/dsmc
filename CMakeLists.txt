cmake_minimum_required(VERSION 3.24)
project(dsmc)

set(CMAKE_CXX_STANDARD 17)

add_executable(dsmc main.cpp Cell.cpp Cell.h Particle.h Grid.cpp Grid.h)

# Mac:
# cd Release; cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-13 ..; make
