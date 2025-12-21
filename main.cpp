#include <iostream>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <fstream>

#include "Grid.h"
#include <json.hpp>

using json = nlohmann::json;

int main() {

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
//    std::string out_dir = R"(C:\Users\windowsuser\iCloudDrive\Projects\coding\monte-carlo\dsmc\OUT\)";
    std::string out_dir = R"(C:\Users\windowsuser\Documents\DSMC_OUT\)";
#else
    // std::string out_dir = "/Users/ian/Downloads/DSMC_OUT/";
    std::string out_dir = "../DSMC_OUT/";
    std::string cmd = "rm " + out_dir + "*.csv";
    std::system(cmd.c_str());
#endif
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

    Grid grid;

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

    grid.d = diam; // kinetic particle diameter
    grid.N = num_particles;
    grid.N_ef = (density/mass)*pow(L,3)/grid.N;
    std::cout << "Each particle represents " << grid.N_ef << " molecules/atoms\n";
    grid.num_dens = density;
    grid.V = L*L*L;
    std::cout << "System volume (V): " << grid.V << std::endl;

    grid.cell_length = L / ncell;

    double v_init = sqrt(3*boltz*T/mass);
    grid.v_mult = v_init;

    grid.a = a_timestep;
    grid.v_max = v_max_mult * v_init; // max particle velocity
    grid.num_dt = config["simulation"]["num_timesteps"];

    grid.create();
    grid.writeParticlesToDisk(out_dir + std::string("afterCreate-") + out_file);

    std::cout << "Particle width (d): " << grid.d << "\n";
    std::cout << "Mean free path: " << grid.lambda << std::endl;
    std::cout << "Cell length: " << grid.cell_length << "\n";
    std::cout << "Grid dimension (for cube): " << grid.dim << "\n";
    std::cout << "Characteristic length: " << grid.charlen<< std::endl;
    std::cout << "Average speed: " << grid.mean_v << std::endl;
    std::cout << "Timestep: " << grid.delta_t << std::endl;


    int write_iter = 0;
    for (int iter=0; iter < grid.num_dt; iter++) {
        grid.calculateCollisionsRejectionSampling();
        grid.updatePositions();
        grid.enforceDomain();
        grid.reassignParticlesToCells();
        std::string iter_str = std::to_string(write_iter);
        auto new_str = std::string(n_zero - std::min(n_zero, (int)iter_str.length()), '0') + iter_str;
        if (iter % skip_every == 0) {
            grid.writeParticlesToDisk(out_dir + new_str + std::string("-") + out_file);
            write_iter += 1;
        }
    }

    return 0;
}
