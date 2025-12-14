//
// MetalCompute.mm
// Metal compute wrapper implementation for DSMC
//

#import "MetalCompute.h"
#import <Foundation/Foundation.h>
#include <iostream>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <cassert>

MetalCompute::MetalCompute()
    : m_device(nil)
    , m_commandQueue(nil)
    , m_library(nil)
    , m_numParticles(0)
    , m_numCells(0)
    , m_hasMesh(false)
{
}

MetalCompute::~MetalCompute() {
    // ARC handles memory management
}

bool MetalCompute::initialize() {
    // Get default Metal device
    m_device = MTLCreateSystemDefaultDevice();
    if (!m_device) {
        std::cerr << "Metal is not supported on this device" << std::endl;
        return false;
    }
    
    std::cout << "Using Metal device: " << [m_device.name UTF8String] << std::endl;
    
    // Create command queue
    m_commandQueue = [m_device newCommandQueue];
    if (!m_commandQueue) {
        std::cerr << "Failed to create command queue" << std::endl;
        return false;
    }
    
    return true;
}

bool MetalCompute::loadShaders(const std::string& metalLibPath) {
    NSError* error = nil;
    
    // Try to load from compiled metallib first
    NSString* libPath = [NSString stringWithUTF8String:metalLibPath.c_str()];
    NSURL* libURL = [NSURL fileURLWithPath:libPath];
    
    m_library = [m_device newLibraryWithURL:libURL error:&error];
    
    if (!m_library) {
        // Fall back to compiling from source
        std::cout << "Metallib not found, compiling from source..." << std::endl;
        
        // Try to find and compile the .metal source file
        NSString* sourcePath = [[libPath stringByDeletingLastPathComponent]
                                stringByAppendingPathComponent:@"Shaders.metal"];
        
        // Also try Metal subdirectory
        if (![[NSFileManager defaultManager] fileExistsAtPath:sourcePath]) {
            sourcePath = [[[libPath stringByDeletingLastPathComponent]
                          stringByAppendingPathComponent:@"Metal"]
                          stringByAppendingPathComponent:@"Shaders.metal"];
        }
        
        NSString* source = [NSString stringWithContentsOfFile:sourcePath
                                                    encoding:NSUTF8StringEncoding
                                                       error:&error];
        if (!source) {
            std::cerr << "Failed to load shader source: " << [[error localizedDescription] UTF8String] << std::endl;
            return false;
        }
        
        MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
        options.mathMode = MTLMathModeFast;
        
        m_library = [m_device newLibraryWithSource:source options:options error:&error];
        if (!m_library) {
            std::cerr << "Failed to compile shaders: " << [[error localizedDescription] UTF8String] << std::endl;
            return false;
        }
    }
    
    std::cout << "Shaders loaded successfully" << std::endl;
    
    // Create compute pipelines
    m_updatePositionsPipeline = createPipeline("updatePositions");
    m_enforceDomainPipeline = createPipeline("enforceDomain");
    m_computeCellIndicesPipeline = createPipeline("computeCellIndices");
    m_countParticlesPerCellPipeline = createPipeline("countParticlesPerCell");
    m_computeCellOffsetsPipeline = createPipeline("computeCellOffsets");
    m_reorderParticlesPipeline = createPipeline("reorderParticles");
    m_calculateCollisionsPipeline = createPipeline("calculateCollisions");
    m_intersectMeshPipeline = createPipeline("intersectMeshPrimitive");
    m_initializeParticlesPipeline = createPipeline("initializeParticles");
    m_copyParticlePositionsPipeline = createPipeline("copyParticlePositions");
    m_clearBufferPipeline = createPipeline("clearBuffer");
    
    return m_updatePositionsPipeline && m_enforceDomainPipeline;
}

id<MTLComputePipelineState> MetalCompute::createPipeline(const std::string& functionName) {
    NSError* error = nil;
    
    NSString* funcName = [NSString stringWithUTF8String:functionName.c_str()];
    id<MTLFunction> function = [m_library newFunctionWithName:funcName];
    
    if (!function) {
        std::cerr << "Failed to find function: " << functionName << std::endl;
        return nil;
    }
    
    id<MTLComputePipelineState> pipeline = [m_device newComputePipelineStateWithFunction:function error:&error];
    
    if (!pipeline) {
        std::cerr << "Failed to create pipeline for " << functionName << ": "
                  << [[error localizedDescription] UTF8String] << std::endl;
        return nil;
    }
    
    return pipeline;
}

void MetalCompute::createParticleBuffers(uint32_t numParticles) {
    m_numParticles = numParticles;
    
    size_t particleBufferSize = numParticles * sizeof(GPUParticle);
    size_t positionsBufferSize = numParticles * sizeof(float) * 4;  // float4 per particle
    size_t indexBufferSize = numParticles * sizeof(uint32_t);
    
    m_particleBuffer = [m_device newBufferWithLength:particleBufferSize
                                             options:MTLResourceStorageModeShared];
    m_particleBufferSorted = [m_device newBufferWithLength:particleBufferSize
                                                   options:MTLResourceStorageModeShared];
    m_positionsBuffer = [m_device newBufferWithLength:positionsBufferSize
                                              options:MTLResourceStorageModeShared];
    m_cellIndicesBuffer = [m_device newBufferWithLength:indexBufferSize
                                                options:MTLResourceStorageModeShared];
    m_particleIndicesBuffer = [m_device newBufferWithLength:indexBufferSize
                                                    options:MTLResourceStorageModeShared];
    
    std::cout << "Created particle buffers for " << numParticles << " particles" << std::endl;
}

void MetalCompute::createCellBuffers(uint32_t numCells) {
    m_numCells = numCells;
    
    size_t cellBufferSize = numCells * sizeof(uint32_t);
    
    m_cellCountsBuffer = [m_device newBufferWithLength:cellBufferSize
                                               options:MTLResourceStorageModeShared];
    m_cellStartsBuffer = [m_device newBufferWithLength:cellBufferSize
                                               options:MTLResourceStorageModeShared];
    m_cellEndsBuffer = [m_device newBufferWithLength:cellBufferSize
                                             options:MTLResourceStorageModeShared];
    m_cellWriteOffsetsBuffer = [m_device newBufferWithLength:cellBufferSize
                                                     options:MTLResourceStorageModeShared];
    m_prefixSumBuffer = [m_device newBufferWithLength:cellBufferSize
                                              options:MTLResourceStorageModeShared];
    
    std::cout << "Created cell buffers for " << numCells << " cells" << std::endl;
}

void MetalCompute::createMeshBuffers(const std::vector<Triangle>& triangles) {
    if (triangles.empty()) {
        m_hasMesh = false;
        return;
    }
    
    size_t triangleBufferSize = triangles.size() * sizeof(Triangle);
    size_t normalBufferSize = triangles.size() * sizeof(float) * 3;
    
    m_triangleBuffer = [m_device newBufferWithLength:triangleBufferSize
                                             options:MTLResourceStorageModeShared];
    m_triangleNormalsBuffer = [m_device newBufferWithLength:normalBufferSize
                                                    options:MTLResourceStorageModeShared];
    
    // Copy triangle data
    memcpy(m_triangleBuffer.contents, triangles.data(), triangleBufferSize);
    
    // Extract and copy normals
    float* normals = (float*)m_triangleNormalsBuffer.contents;
    for (size_t i = 0; i < triangles.size(); i++) {
        normals[i * 3 + 0] = triangles[i].normal[0];
        normals[i * 3 + 1] = triangles[i].normal[1];
        normals[i * 3 + 2] = triangles[i].normal[2];
    }
    
    std::cout << "Created mesh buffers for " << triangles.size() << " triangles" << std::endl;
}

bool MetalCompute::buildAccelerationStructure(const std::vector<Triangle>& triangles) {
    if (triangles.empty()) {
        m_hasMesh = false;
        return true;
    }
    
    // Create geometry descriptor for triangles
    MTLAccelerationStructureTriangleGeometryDescriptor* geometryDesc =
        [MTLAccelerationStructureTriangleGeometryDescriptor descriptor];
    
    // Create vertex buffer (3 vertices per triangle, 3 floats per vertex)
    size_t vertexBufferSize = triangles.size() * 3 * 3 * sizeof(float);
    id<MTLBuffer> vertexBuffer = [m_device newBufferWithLength:vertexBufferSize
                                                       options:MTLResourceStorageModeShared];
    
    float* vertices = (float*)vertexBuffer.contents;
    for (size_t i = 0; i < triangles.size(); i++) {
        // Vertex 0
        vertices[i * 9 + 0] = triangles[i].v0[0];
        vertices[i * 9 + 1] = triangles[i].v0[1];
        vertices[i * 9 + 2] = triangles[i].v0[2];
        // Vertex 1
        vertices[i * 9 + 3] = triangles[i].v1[0];
        vertices[i * 9 + 4] = triangles[i].v1[1];
        vertices[i * 9 + 5] = triangles[i].v1[2];
        // Vertex 2
        vertices[i * 9 + 6] = triangles[i].v2[0];
        vertices[i * 9 + 7] = triangles[i].v2[1];
        vertices[i * 9 + 8] = triangles[i].v2[2];
    }
    
    geometryDesc.vertexBuffer = vertexBuffer;
    geometryDesc.vertexBufferOffset = 0;
    geometryDesc.vertexStride = sizeof(float) * 3;
    geometryDesc.triangleCount = triangles.size();
    geometryDesc.indexBuffer = nil;  // Non-indexed triangles
    
    // Create primitive acceleration structure descriptor
    MTLPrimitiveAccelerationStructureDescriptor* accelDesc =
        [MTLPrimitiveAccelerationStructureDescriptor descriptor];
    accelDesc.geometryDescriptors = @[geometryDesc];
    
    // Get sizes needed for acceleration structure
    MTLAccelerationStructureSizes sizes = [m_device accelerationStructureSizesWithDescriptor:accelDesc];
    
    // Create acceleration structure
    m_accelerationStructure = [m_device newAccelerationStructureWithSize:sizes.accelerationStructureSize];
    
    // Create scratch buffer for building
    id<MTLBuffer> scratchBuffer = [m_device newBufferWithLength:sizes.buildScratchBufferSize
                                                        options:MTLResourceStorageModePrivate];
    
    // Build acceleration structure
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLAccelerationStructureCommandEncoder> encoder = [commandBuffer accelerationStructureCommandEncoder];
    
    [encoder buildAccelerationStructure:m_accelerationStructure
                             descriptor:accelDesc
                          scratchBuffer:scratchBuffer
                    scratchBufferOffset:0];
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    
    m_hasMesh = true;
    std::cout << "Built acceleration structure for " << triangles.size() << " triangles" << std::endl;
    
    return true;
}

void MetalCompute::uploadParticles(const std::vector<GPUParticle>& particles) {
    if (particles.size() != m_numParticles) {
        std::cerr << "Particle count mismatch!" << std::endl;
        return;
    }
    
    memcpy(m_particleBuffer.contents, particles.data(), particles.size() * sizeof(GPUParticle));
}

void MetalCompute::downloadParticles(std::vector<GPUParticle>& particles) {
    particles.resize(m_numParticles);
    memcpy(particles.data(), m_particleBuffer.contents, m_numParticles * sizeof(GPUParticle));
}

void MetalCompute::downloadPositions(std::vector<float>& positions) {
    positions.resize(m_numParticles * 4);
    
    // Run kernel to copy positions
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_copyParticlePositionsPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_positionsBuffer offset:0 atIndex:1];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:2];
    
    dispatchCompute(encoder, m_copyParticlePositionsPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    
    memcpy(positions.data(), m_positionsBuffer.contents, m_numParticles * 4 * sizeof(float));
}

void MetalCompute::setSimulationParams(const SimulationParams& params) {
    m_params = params;
    
    if (!m_paramsBuffer) {
        m_paramsBuffer = [m_device newBufferWithLength:sizeof(SimulationParams)
                                               options:MTLResourceStorageModeShared];
    }
    
    memcpy(m_paramsBuffer.contents, &m_params, sizeof(SimulationParams));
}

void MetalCompute::dispatchCompute(id<MTLComputeCommandEncoder> encoder,
                                   id<MTLComputePipelineState> pipeline,
                                   uint32_t numThreads) {
    NSUInteger threadGroupSize = pipeline.maxTotalThreadsPerThreadgroup;
    if (threadGroupSize > 256) threadGroupSize = 256;
    
    MTLSize threadsPerGroup = MTLSizeMake(threadGroupSize, 1, 1);
    MTLSize numGroups = MTLSizeMake((numThreads + threadGroupSize - 1) / threadGroupSize, 1, 1);
    
    [encoder dispatchThreadgroups:numGroups threadsPerThreadgroup:threadsPerGroup];
}

void MetalCompute::initializeParticles(uint32_t seed) {
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_initializeParticlesPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:1];
    [encoder setBytes:&seed length:sizeof(uint32_t) atIndex:2];
    
    dispatchCompute(encoder, m_initializeParticlesPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalCompute::updatePositions() {
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_updatePositionsPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:1];
    
    dispatchCompute(encoder, m_updatePositionsPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalCompute::enforceDomain() {
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_enforceDomainPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:1];
    
    dispatchCompute(encoder, m_enforceDomainPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalCompute::computeCellIndices() {
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_computeCellIndicesPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_cellIndicesBuffer offset:0 atIndex:1];
    [encoder setBuffer:m_particleIndicesBuffer offset:0 atIndex:2];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:3];
    
    dispatchCompute(encoder, m_computeCellIndicesPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalCompute::countParticlesPerCell() {
    // Clear cell counts first
    id<MTLCommandBuffer> clearBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> clearEncoder = [clearBuffer computeCommandEncoder];
    
    [clearEncoder setComputePipelineState:m_clearBufferPipeline];
    [clearEncoder setBuffer:m_cellCountsBuffer offset:0 atIndex:0];
    [clearEncoder setBytes:&m_numCells length:sizeof(uint32_t) atIndex:1];
    
    dispatchCompute(clearEncoder, m_clearBufferPipeline, m_numCells);
    
    [clearEncoder endEncoding];
    [clearBuffer commit];
    [clearBuffer waitUntilCompleted];
    
    // Count particles per cell
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_countParticlesPerCellPipeline];
    [encoder setBuffer:m_cellIndicesBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_cellCountsBuffer offset:0 atIndex:1];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:2];
    
    dispatchCompute(encoder, m_countParticlesPerCellPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalCompute::computePrefixSum() {
    // Simple CPU-side prefix sum (exclusive scan)
    uint32_t* counts = (uint32_t*)m_cellCountsBuffer.contents;
    uint32_t* prefixSum = (uint32_t*)m_prefixSumBuffer.contents;
    uint32_t* starts = (uint32_t*)m_cellStartsBuffer.contents;
    uint32_t* ends = (uint32_t*)m_cellEndsBuffer.contents;
    
    uint32_t sum = 0;
    for (uint32_t i = 0; i < m_numCells; i++) {
        starts[i] = sum;
        sum += counts[i];
        ends[i] = sum;
        prefixSum[i] = sum;
    }
}

void MetalCompute::reorderParticles() {
    // Reorder particles (Gather)
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_reorderParticlesPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_particleBufferSorted offset:0 atIndex:1];
    [encoder setBuffer:m_particleIndicesBuffer offset:0 atIndex:2];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:3];
    
    dispatchCompute(encoder, m_reorderParticlesPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    
    // Swap buffers
    id<MTLBuffer> temp = m_particleBuffer;
    m_particleBuffer = m_particleBufferSorted;
    m_particleBufferSorted = temp;
}

void MetalCompute::calculateCollisions(uint32_t frameNumber) {
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_calculateCollisionsPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_cellStartsBuffer offset:0 atIndex:1];
    [encoder setBuffer:m_cellEndsBuffer offset:0 atIndex:2];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:3];
    [encoder setBytes:&frameNumber length:sizeof(uint32_t) atIndex:4];
    
    // One thread per cell for collision processing
    dispatchCompute(encoder, m_calculateCollisionsPipeline, m_numCells);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalCompute::intersectMesh(uint32_t frameNumber) {
    if (!m_hasMesh || !m_accelerationStructure) return;
    
    id<MTLCommandBuffer> commandBuffer = [m_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:m_intersectMeshPipeline];
    [encoder setBuffer:m_particleBuffer offset:0 atIndex:0];
    [encoder setBuffer:m_paramsBuffer offset:0 atIndex:1];
    [encoder setAccelerationStructure:m_accelerationStructure atBufferIndex:2];
    [encoder setBuffer:m_triangleNormalsBuffer offset:0 atIndex:3];
    [encoder setBytes:&frameNumber length:sizeof(uint32_t) atIndex:4];
    
    dispatchCompute(encoder, m_intersectMeshPipeline, m_numParticles);
    
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalCompute::runSimulationStep(uint32_t frameNumber, bool hasMesh) {
    // Check initial state
    checkForNaNs(frameNumber, "start_of_frame");
    
    // Phase 1: Collisions (must happen before position update for proper DSMC)
    computeCellIndices();
    
    // CPU Sort
    {
        // Download cell indices
        uint32_t* cellIndicesGPU = (uint32_t*)m_cellIndicesBuffer.contents;
        std::vector<uint32_t> cellIndices(m_numParticles);
        memcpy(cellIndices.data(), cellIndicesGPU, m_numParticles * sizeof(uint32_t));
        
        // Create and sort particle indices
        std::vector<uint32_t> particleIndices(m_numParticles);
        std::iota(particleIndices.begin(), particleIndices.end(), 0);
        
        std::stable_sort(particleIndices.begin(), particleIndices.end(),
                         [&](uint32_t a, uint32_t b) {
                             return cellIndices[a] < cellIndices[b];
                         });
                         
        // Upload sorted particle indices
        memcpy(m_particleIndicesBuffer.contents, particleIndices.data(), m_numParticles * sizeof(uint32_t));
        
        // Compute cell starts/ends on CPU
        uint32_t* starts = (uint32_t*)m_cellStartsBuffer.contents;
        uint32_t* ends = (uint32_t*)m_cellEndsBuffer.contents;
        
        // Initialize to 0
        memset(starts, 0, m_numCells * sizeof(uint32_t));
        memset(ends, 0, m_numCells * sizeof(uint32_t));
        
        if (m_numParticles > 0) {
            uint32_t currentCell = cellIndices[particleIndices[0]];
            starts[currentCell] = 0;
            
            for (uint32_t i = 1; i < m_numParticles; i++) {
                uint32_t cell = cellIndices[particleIndices[i]];
                if (cell != currentCell) {
                    ends[currentCell] = i;
                    starts[cell] = i;
                    currentCell = cell;
                }
            }
            ends[currentCell] = m_numParticles;
        }
    }
    
    // Reorder (Gather)
    reorderParticles();
    checkForNaNs(frameNumber, "after_reorderParticles");
    
    calculateCollisions(frameNumber);
    checkForNaNs(frameNumber, "calculateCollisions");
    
    // Phase 2: Position update
    updatePositions();
    checkForNaNs(frameNumber, "updatePositions");
    
    // Phase 3: Mesh intersection (if mesh exists)
    if (hasMesh && m_hasMesh) {
        intersectMesh(frameNumber);
        checkForNaNs(frameNumber, "intersectMesh");
    }
    
    // Phase 4: Domain boundary enforcement
    enforceDomain();
    checkForNaNs(frameNumber, "enforceDomain");
}

void MetalCompute::checkForNaNs(uint32_t frameNumber, const char* phase) {
    GPUParticle* particles = (GPUParticle*)m_particleBuffer.contents;
    
    for (uint32_t i = 0; i < m_numParticles; i++) {
        // Check positions
        if (std::isnan(particles[i].pos[0]) ||
            std::isnan(particles[i].pos[1]) ||
            std::isnan(particles[i].pos[2])) {
            std::cerr << "Particle " << i << " pos: ("
                      << particles[i].pos[0] << ", "
                      << particles[i].pos[1] << ", "
                      << particles[i].pos[2] << ")" << std::endl;
            std::cerr << "Particle " << i << " vel: ("
                      << particles[i].vel[0] << ", "
                      << particles[i].vel[1] << ", "
                      << particles[i].vel[2] << ")" << std::endl;
            throw std::runtime_error("NaN detected in particle " + std::to_string(i) +
                                   " position at frame " + std::to_string(frameNumber) +
                                   " after " + phase);
        }
        
        // Check velocities
        if (std::isnan(particles[i].vel[0]) ||
            std::isnan(particles[i].vel[1]) ||
            std::isnan(particles[i].vel[2])) {
            std::cerr << "Particle " << i << " pos: ("
                      << particles[i].pos[0] << ", "
                      << particles[i].pos[1] << ", "
                      << particles[i].pos[2] << ")" << std::endl;
            std::cerr << "Particle " << i << " vel: ("
                      << particles[i].vel[0] << ", "
                      << particles[i].vel[1] << ", "
                      << particles[i].vel[2] << ")" << std::endl;
            throw std::runtime_error("NaN detected in particle " + std::to_string(i) +
                                   " velocity at frame " + std::to_string(frameNumber) +
                                   " after " + phase);
        }
    }
}

void MetalCompute::waitForCompletion() {
    // All operations already wait for completion due to waitUntilCompleted calls
}
