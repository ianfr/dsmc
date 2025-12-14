//
// Grid.mm
// Metal-accelerated DSMC Grid implementation
//

#include "Grid.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

#ifndef PI
#define PI 3.14159265
#endif

Grid::Grid()
    : N(0)
    , N_ef(1)
    , V(0)
    , num_dens(0)
    , d(0)
    , cell_length(0)
    , v_mult(0)
    , v_max(0)
    , a(0)
    , num_dt(0)
    , lambda(0)
    , charlen(0)
    , mean_v(0)
    , delta_t(0)
    , dim(0)
    , m_initialized(false)
    , m_hasMesh(false)
{
    m_metal = std::make_unique<MetalCompute>();
}

Grid::~Grid() = default;

bool Grid::initializeMetal(const std::string& shaderPath) {
    if (!m_metal->initialize()) {
        std::cerr << "Failed to initialize Metal" << std::endl;
        return false;
    }
    
    if (!m_metal->loadShaders(shaderPath)) {
        std::cerr << "Failed to load Metal shaders" << std::endl;
        return false;
    }
    
    m_initialized = true;
    return true;
}

void Grid::setupSimulationParams() {
    SimulationParams params;
    params.delta_t = static_cast<float>(delta_t);
    params.v_max = static_cast<float>(v_max);
    params.d = static_cast<float>(d);
    params.cell_length = static_cast<float>(cell_length);
    params.domain_min = 0.0f;
    params.domain_max = static_cast<float>(dim * cell_length);
    params.v_mult = static_cast<float>(v_mult);
    params.eps = 1e-8f;
    params.num_particles = static_cast<uint32_t>(N);
    params.num_cells_per_dim = static_cast<uint32_t>(dim);
    params.f_n = static_cast<uint32_t>(N_ef);
    params.V_c = static_cast<float>(cell_length * cell_length * cell_length);
    
    m_metal->setSimulationParams(params);
}

void Grid::create() {
    if (!m_initialized) {
        std::cerr << "Metal not initialized! Call initializeMetal() first." << std::endl;
        return;
    }
    
    // Compute derived values
    lambda = V / (sqrt(2) * N_ef * PI * d * d);
    charlen = std::pow((cell_length * V) / N_ef, 1.0 / 3.0);
    dim = static_cast<int>(pow(V, 1.0 / 3.0) / cell_length);
    
    // Estimate mean velocity for timestep calculation
    mean_v = v_mult;  // Initial estimate
    delta_t = (a * cell_length) / mean_v;
    
    std::cout << "Grid configuration:" << std::endl;
    std::cout << "  Particles: " << N << std::endl;
    std::cout << "  Particles per simulated: " << N_ef << std::endl;
    std::cout << "  Grid dimension: " << dim << "x" << dim << "x" << dim << std::endl;
    std::cout << "  Cell length: " << cell_length << std::endl;
    std::cout << "  Domain size: " << dim * cell_length << std::endl;
    std::cout << "  Mean free path: " << lambda << std::endl;
    std::cout << "  Timestep: " << delta_t << std::endl;
    
    // Create GPU buffers
    m_metal->createParticleBuffers(static_cast<uint32_t>(N));
    m_metal->createCellBuffers(static_cast<uint32_t>(dim * dim * dim));
    
    // Setup simulation parameters
    setupSimulationParams();
    
    // Initialize particles on GPU with random positions and velocities
    std::random_device rd;
    m_metal->initializeParticles(rd());
    
    std::cout << "Grid created and particles initialized on GPU" << std::endl;
}

void Grid::loadMesh(const std::string& meshPath) {
    Mesh mesh;
    if (!mesh.load(meshPath)) {
        std::cerr << "Failed to load mesh: " << meshPath << std::endl;
        return;
    }
    
    loadMesh(mesh);
}

void Grid::loadMesh(const Mesh& mesh) {
    if (mesh.getTriangleCount() == 0) {
        std::cerr << "Mesh has no triangles" << std::endl;
        return;
    }
    
    // Fit mesh to domain (with some margin)
    Mesh transformedMesh = mesh;
    float domainSize = static_cast<float>(dim * cell_length);
    float margin = domainSize * 0.1f;  // 10% margin
    transformedMesh.fitToBounds(
        margin, domainSize - margin,
        margin, domainSize - margin,
        margin, domainSize - margin
    );
    
    // Create mesh buffers and acceleration structure
    m_metal->createMeshBuffers(transformedMesh.getTriangles());
    
    if (m_metal->buildAccelerationStructure(transformedMesh.getTriangles())) {
        m_hasMesh = true;
        std::cout << "Mesh loaded with " << transformedMesh.getTriangleCount()
                  << " triangles" << std::endl;
    } else {
        std::cerr << "Failed to build acceleration structure" << std::endl;
    }
}

void Grid::createTestSphere(float radius) {
    float domainSize = static_cast<float>(dim * cell_length);
    float center = domainSize / 2.0f;
    
    Mesh sphere = Mesh::createSphere(radius, 3);  // 3 subdivisions = 1280 triangles
    sphere.translate(center, center, center);
    
    m_metal->createMeshBuffers(sphere.getTriangles());
    
    if (m_metal->buildAccelerationStructure(sphere.getTriangles())) {
        m_hasMesh = true;
        std::cout << "Test sphere created with " << sphere.getTriangleCount()
                  << " triangles at center (" << center << ", " << center << ", " << center << ")"
                  << std::endl;
    }
}

void Grid::createTestBox(float width, float height, float depth) {
    float domainSize = static_cast<float>(dim * cell_length);
    float center = domainSize / 2.0f;
    
    Mesh box = Mesh::createBox(width, height, depth);
    box.translate(center, center, center);
    
    m_metal->createMeshBuffers(box.getTriangles());
    
    if (m_metal->buildAccelerationStructure(box.getTriangles())) {
        m_hasMesh = true;
        std::cout << "Test box created at center" << std::endl;
    }
}

void Grid::runSimulationStep(uint32_t frameNumber) {
    if (!m_initialized) {
        std::cerr << "Metal not initialized!" << std::endl;
        return;
    }
    
    m_metal->runSimulationStep(frameNumber, m_hasMesh);
}

void Grid::getParticlePositions(std::vector<float>& positions) {
    m_metal->downloadPositions(positions);
}

void Grid::writeParticlesToDisk(const std::string& filename) {
    // Download particle positions from GPU
    std::vector<float> positions;
    m_metal->downloadPositions(positions);
    
    // Write to file
    std::ofstream outFile(filename, std::ios::trunc);
    outFile << "x, y, z, velmag, cell\n";
    
    uint32_t numParticles = m_metal->getNumParticles();
    float cellLen = static_cast<float>(cell_length);
    uint32_t numCellsPerDim = static_cast<uint32_t>(dim);
    
    for (uint32_t i = 0; i < numParticles; i++) {
        float x = positions[i * 4 + 0];
        float y = positions[i * 4 + 1];
        float z = positions[i * 4 + 2];
        float velmag = positions[i * 4 + 3];
        
        // Compute cell index from current position
        uint32_t ix = std::min(static_cast<uint32_t>(x / cellLen), numCellsPerDim - 1);
        uint32_t iy = std::min(static_cast<uint32_t>(y / cellLen), numCellsPerDim - 1);
        uint32_t iz = std::min(static_cast<uint32_t>(z / cellLen), numCellsPerDim - 1);
        uint32_t cellIdx = ix * numCellsPerDim * numCellsPerDim + iy * numCellsPerDim + iz;
        
        outFile << x << ", " << y << ", " << z << ", " << velmag << ", " << cellIdx << "\n";
    }
    
    outFile.close();
}
