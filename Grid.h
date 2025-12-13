//
// Grid.h
// Metal-accelerated DSMC Grid
//

#ifndef DSMC_GRID_H
#define DSMC_GRID_H

#include "MetalCompute.h"
#include "Mesh.h"
#include <string>
#include <vector>
#include <memory>

class Grid {
public:
    Grid();
    ~Grid();
    
    // Configuration (set before calling create())
    int N;                  // number of simulated particles
    int N_ef;               // particles represented by each simulated particle
    double V;               // system volume
    double num_dens;        // number density
    double d;               // particle diameter
    double cell_length;     // cell side length
    double v_mult;          // initial velocity multiplier
    double v_max;           // maximum velocity
    double a;               // timestep multiplier
    int num_dt;             // number of timesteps
    
    // Computed values
    double lambda;          // mean free path
    double charlen;         // characteristic length
    double mean_v;          // mean velocity
    double delta_t;         // timestep
    int dim;                // grid dimension (cells per side)
    
    // Methods
    bool initializeMetal(const std::string& shaderPath);
    void create();
    void loadMesh(const std::string& meshPath);
    void loadMesh(const Mesh& mesh);
    void createTestSphere(float radius);
    void createTestBox(float width, float height, float depth);
    
    // Simulation step (all on GPU)
    void runSimulationStep(uint32_t frameNumber);
    
    // Output
    void writeParticlesToDisk(const std::string& filename);
    
    // Get particle data for visualization
    void getParticlePositions(std::vector<float>& positions);
    
private:
    std::unique_ptr<MetalCompute> m_metal;
    bool m_initialized;
    bool m_hasMesh;
    
    // CPU-side particle storage for output
    std::vector<GPUParticle> m_particles;
    
    void setupSimulationParams();
};

#endif // DSMC_GRID_H
