cmake_minimum_required(VERSION 3.24)
project(dsmc)

set(CMAKE_CXX_STANDARD 17)

enable_language(OBJCXX)
find_library(METAL_LIBRARY Metal REQUIRED)
find_library(METALKIT_LIBRARY MetalKit REQUIRED)
find_library(FOUNDATION_LIBRARY Foundation REQUIRED)

# Compile Metal shaders
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/Shaders.metallib
    COMMAND xcrun -sdk macosx metal -c ${CMAKE_SOURCE_DIR}/Metal/Shaders.metal -o ${CMAKE_BINARY_DIR}/Shaders.air
    COMMAND xcrun -sdk macosx metallib ${CMAKE_BINARY_DIR}/Shaders.air -o ${CMAKE_BINARY_DIR}/Shaders.metallib
    DEPENDS ${CMAKE_SOURCE_DIR}/Metal/Shaders.metal
    COMMENT "Compiling Metal shaders"
)
add_custom_target(MetalShaders ALL DEPENDS ${CMAKE_BINARY_DIR}/Shaders.metallib)

include_directories(
    ${CMAKE_SOURCE_DIR}
)

# Metal-accelerated executable
set_source_files_properties(
    main.cpp 
    MetalCompute.mm MetalCompute.h 
    Grid.mm Grid.h 
    Mesh.cpp Mesh.h
    PROPERTIES LANGUAGE OBJCXX
)
add_executable(dsmc 
    main.cpp 
    MetalCompute.mm MetalCompute.h 
    Grid.mm Grid.h 
    Mesh.cpp Mesh.h
)
target_link_libraries(dsmc ${METAL_LIBRARY} ${METALKIT_LIBRARY} ${FOUNDATION_LIBRARY})
add_dependencies(dsmc MetalShaders)

# Copy shader library to build directory
add_custom_command(TARGET dsmc POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/Shaders.metallib
        $<TARGET_FILE_DIR:dsmc>/Shaders.metallib
)
