//
// main.cpp
// Metal-accelerated DSMC simulation
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <chrono>

#include "Grid.h"
#include "json.hpp"

using json = nlohmann::json;

int main(int argc, char* argv[]) {
    std::cout << "DSMC - Metal Accelerated Version" << std::endl;
    std::cout << "=================================" << std::endl;

    std::string out_dir = "../DSMC_OUT/";
    std::string cmd = "rm -f " + out_dir + "*.csv";
    std::system(cmd.c_str());

    std::string out_file = "particle.csv";

    // Load configuration from JSON file
    std::ifstream config_file("../config.json");
    if (!config_file.is_open()) {
        std::cerr << "Error: Could not open config.json" << std::endl;
        return 1;
    }
    json config = json::parse(config_file);
    config_file.close();

    // Read simulation parameters
    int skip_every = config["simulation"]["skip_every"];
    int n_zero = config["simulation"]["n_zero"];

    // Read physical constants and gas properties
    double boltz = config["physical_constants"]["boltzmann"];
    double mass = config["gas_properties"]["mass"];
    double diam = config["gas_properties"]["diameter"];
    double T = config["gas_properties"]["temperature"];
    double density = config["gas_properties"]["density"];
    double L = config["domain"]["length"];

    // Read particle and domain settings
    int ncell = config["domain"]["num_cells_per_dim"];
    int num_particles = config["particles"]["num_particles"];
    double a_timestep = config["particles"]["a_timestep"];
    double v_max_mult = config["particles"]["v_max_multiplier"];

    // Optional mesh configuration
    std::string mesh_file = "";
    std::string mesh_type = "none";
    float mesh_size = 0.0f;
    
    if (config.contains("mesh")) {
        if (config["mesh"].contains("file")) {
            mesh_file = config["mesh"]["file"];
        }
        if (config["mesh"].contains("type")) {
            mesh_type = config["mesh"]["type"];
        }
        if (config["mesh"].contains("size")) {
            mesh_size = config["mesh"]["size"];
        }
    }

    // Create Metal-accelerated grid
    Grid grid;

    // Set grid parameters
    grid.d = diam;
    grid.N = num_particles;
    grid.N_ef = static_cast<int>((density / mass) * pow(L, 3) / grid.N);
    std::cout << "Each particle represents " << grid.N_ef << " molecules/atoms" << std::endl;
    
    grid.num_dens = density;
    grid.V = L * L * L;
    std::cout << "System volume (V): " << grid.V << std::endl;

    grid.cell_length = L / ncell;

    double v_init = sqrt(3 * boltz * T / mass);
    grid.v_mult = v_init;

    grid.a = a_timestep;
    grid.v_max = v_max_mult * v_init;
    grid.num_dt = config["simulation"]["num_timesteps"];

    // Initialize Metal
    std::string shader_path = "./Shaders.metallib";
    
    if (!grid.initializeMetal(shader_path)) {
        std::cerr << "Failed to initialize Metal. Exiting." << std::endl;
        return 1;
    }

    // Create grid and particles
    grid.create();

    std::cout << "\nSimulation parameters:" << std::endl;
    std::cout << "  Particle diameter: " << grid.d << std::endl;
    std::cout << "  Mean free path: " << grid.lambda << std::endl;
    std::cout << "  Cell length: " << grid.cell_length << std::endl;
    std::cout << "  Grid dimension: " << grid.dim << std::endl;
    std::cout << "  Characteristic length: " << grid.charlen << std::endl;
    std::cout << "  Average speed: " << grid.mean_v << std::endl;
    std::cout << "  Timestep: " << grid.delta_t << std::endl;

    // Load or create mesh if specified
    if (!mesh_file.empty()) {
        std::cout << "\nLoading mesh from: " << mesh_file << std::endl;
        grid.loadMesh(mesh_file);
    } else if (mesh_type == "sphere") {
        float radius = (mesh_size > 0) ? mesh_size : static_cast<float>(L * 0.2);
        std::cout << "\nCreating test sphere with radius: " << radius << std::endl;
        grid.createTestSphere(radius);
    } else if (mesh_type == "box") {
        float size = (mesh_size > 0) ? mesh_size : static_cast<float>(L * 0.3);
        std::cout << "\nCreating test box with size: " << size << std::endl;
        grid.createTestBox(size, size, size);
    }

    // Remove any particles inside the mesh if it exists
    grid.pruneParticlesInsideMesh();

    // Write initial state
    grid.writeParticlesToDisk(out_dir + "afterCreate-" + out_file);

    std::cout << "\nStarting simulation (" << grid.num_dt << " timesteps)..." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    
    int write_iter = 0;
    for (int iter = 0; iter < grid.num_dt; iter++) {
        // Run full simulation step on GPU
        grid.runSimulationStep(static_cast<uint32_t>(iter));

        // Write output if needed
        if (iter % skip_every == 0) {
            std::string iter_str = std::to_string(write_iter);
            auto new_str = std::string(n_zero - std::min(n_zero, static_cast<int>(iter_str.length())), '0') + iter_str;
            grid.writeParticlesToDisk(out_dir + new_str + "-" + out_file);
            write_iter++;
        }

        // Progress report
        if ((iter + 1) % 100 == 0 || iter == grid.num_dt - 1) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
            double progress = 100.0 * (iter + 1) / grid.num_dt;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% "
                      << "(" << iter + 1 << "/" << grid.num_dt << ") "
                      << "Elapsed: " << elapsed / 1000.0 << "s" << std::flush;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    std::cout << "\n\nSimulation complete!" << std::endl;
    std::cout << "Total time: " << total_time / 1000.0 << " seconds" << std::endl;
    std::cout << "Average time per step: " << total_time / static_cast<double>(grid.num_dt) << " ms" << std::endl;

    return 0;
}
