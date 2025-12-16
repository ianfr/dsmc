//
// MetalCompute.h
// Metal compute wrapper for DSMC
//

#ifndef DSMC_METALCOMPUTE_H
#define DSMC_METALCOMPUTE_H

#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>

#include <vector>
#include <string>
#include <memory>

#include "Mesh.h"  // For Triangle struct

// Must match the GPU struct layout exactly (float3 aligned to 16 bytes)
struct GPUParticle {
    float pos[3];
    float pad0;
    float vel[3];
    float pad1;
};

struct SimulationParams {
    float delta_t;
    float v_max;
    float d;
    float cell_length;
    float domain_min;
    float domain_max;
    float v_mult;
    float eps;
    uint32_t num_particles;
    uint32_t num_cells_per_dim;
    uint32_t f_n;
    float V_c;
};

class MetalCompute {
public:
    MetalCompute();
    ~MetalCompute();
    
    // Initialization
    bool initialize();
    bool loadShaders(const std::string& metalLibPath);
    
    // Buffer management
    void createParticleBuffers(uint32_t numParticles);
    void createCellBuffers(uint32_t numCells);
    void createMeshBuffers(const std::vector<Triangle>& triangles);
    
    // Upload/download data
    void uploadParticles(const std::vector<GPUParticle>& particles);
    void downloadParticles(std::vector<GPUParticle>& particles);
    void downloadPositions(std::vector<float>& positions);  // x,y,z,velmag per particle
    void downloadCellIndices(std::vector<uint32_t>& cellIndices);  // cell index per particle
    void setSimulationParams(const SimulationParams& params);
    
    // Acceleration structure for ray tracing
    bool buildAccelerationStructure(const std::vector<Triangle>& triangles);
    
    // Compute operations
    void initializeParticles(uint32_t seed);
    void updatePositions();
    void enforceDomain();
    void computeCellIndices();
    void countParticlesPerCell();
    void reorderParticles();
    void calculateCollisions(uint32_t frameNumber);
    void intersectMesh(uint32_t frameNumber);
    uint32_t pruneParticlesInsideMesh(); // Called from Grid::pruneParticlesInsideMesh
    
    // Full simulation step
    void runSimulationStep(uint32_t frameNumber, bool hasMesh);
    
    // Synchronization
    void waitForCompletion();
    
    // NaN checking - throws exception if NaNs detected
    void checkForNaNs(uint32_t frameNumber, const char* phase);
    
    // Getters
    uint32_t getNumParticles() const { return m_numParticles; }
    uint32_t getNumCells() const { return m_numCells; }
    bool hasMesh() const { return m_hasMesh; }
    
private:
    // Metal objects
    id<MTLDevice> m_device;
    id<MTLCommandQueue> m_commandQueue;
    id<MTLLibrary> m_library;
    
    // Compute pipelines
    id<MTLComputePipelineState> m_updatePositionsPipeline;
    id<MTLComputePipelineState> m_updatePositionsWithMeshPipeline;
    id<MTLComputePipelineState> m_enforceDomainPipeline;
    id<MTLComputePipelineState> m_computeCellIndicesPipeline;
    id<MTLComputePipelineState> m_countParticlesPerCellPipeline;
    id<MTLComputePipelineState> m_computeCellOffsetsPipeline;
    id<MTLComputePipelineState> m_reorderParticlesPipeline;
    id<MTLComputePipelineState> m_calculateCollisionsPipeline;
    id<MTLComputePipelineState> m_intersectMeshPipeline;
    id<MTLComputePipelineState> m_pruneParticlesInsideMeshPipeline;
    id<MTLComputePipelineState> m_ejectParticlesFromMeshPipeline;
    id<MTLComputePipelineState> m_initializeParticlesPipeline;
    id<MTLComputePipelineState> m_copyParticlePositionsPipeline;
    id<MTLComputePipelineState> m_clearBufferPipeline;
    
    // Particle buffers (double-buffered for sorting)
    id<MTLBuffer> m_particleBuffer;
    id<MTLBuffer> m_particleBufferSorted;
    id<MTLBuffer> m_positionsBuffer;  // For output (float4: x,y,z,velmag)
    
    // Spatial hashing buffers
    id<MTLBuffer> m_cellIndicesBuffer;
    id<MTLBuffer> m_particleIndicesBuffer;
    id<MTLBuffer> m_cellCountsBuffer;
    id<MTLBuffer> m_cellStartsBuffer;
    id<MTLBuffer> m_cellEndsBuffer;
    id<MTLBuffer> m_cellWriteOffsetsBuffer;
    id<MTLBuffer> m_prefixSumBuffer;
    
    // Mesh and acceleration structure
    id<MTLBuffer> m_triangleBuffer;
    id<MTLBuffer> m_triangleNormalsBuffer;
    id<MTLAccelerationStructure> m_accelerationStructure;
    bool m_hasMesh;
    
    // Simulation parameters
    id<MTLBuffer> m_paramsBuffer;
    SimulationParams m_params;
    uint32_t m_numParticles;
    uint32_t m_numCells;
    
    // Helper methods
    id<MTLComputePipelineState> createPipeline(const std::string& functionName);
    void dispatchCompute(id<MTLComputeCommandEncoder> encoder,
                        id<MTLComputePipelineState> pipeline,
                        uint32_t numThreads);
    void computePrefixSum();  // CPU-side prefix sum for cell offsets
};

#endif // DSMC_METALCOMPUTE_H
