#include <iostream>
#include <cmath>
#include <math.h>
#include <stdlib.h>

#include "Grid.h"

int main() {

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
//    std::string out_dir = R"(C:\Users\windowsuser\iCloudDrive\Projects\coding\monte-carlo\dsmc\OUT\)";
    std::string out_dir = R"(C:\Users\windowsuser\Documents\DSMC_OUT\)";
#else
    std::string out_dir = "../DSMC_OUT/";
    std::string cmd = "rm " + out_dir + "*.csv";
    std::system(cmd.c_str());
#endif
    std::string out_file = "particle.csv";

    int skip_every = 1; // skip every _ steps when writing out
    int n_zero = 10; // 10 digits for indexing output CSV series

    Grid grid;

    // True constants
    constexpr double boltz = 1.3806e-23; // J/K
    constexpr double mass = 6.63e-26; // mass argon
    constexpr double diam = 3.66e-10; // eff diam argon
    constexpr double T = 273; // temperature (K)
    constexpr double density = 1.78; // density of argon at STP (kg/m^3)
    double pack_density_prop = 0.65; // https://en.wikipedia.org/wiki/Random_close_pack

    // Parameters
    // We assume a cubic volume
    // double approx_volume = 1000; // ~10^3
    int n_sim = 1e5;
    int alpha_inv = 10; // the numeber of mean free paths per cell
    double ncell = 100; // number of cells in ONE dimension
    double packet_radius = 0.5e-4; // 1mm packet
    
    // double L = 100 * 1e-6; // x microns
    // double L = 1;

    // double mean_free_path = volume / (std::sqrt(2) * n_sim * M_PI * diam * diam);
    // double L = 0.1 * mean_free_path;

    // Each simulated particle will actually be simulating a collection OF 'packets' OF argon atoms
    double packet_volume = (4/3) * M_PI * std::pow(packet_radius,3);
    double packet_density = density * pack_density_prop;
    double packet_mass = density * packet_volume;
    double packet_diam = packet_radius * 2;

    // double L = std::sqrt(alpha_inv * M_SQRT2 * n_sim * M_PI * diam * diam);
    double L = std::sqrt(alpha_inv * M_SQRT2 * n_sim * M_PI * packet_diam * packet_diam);

    grid.d = packet_diam; // packet kinetic particle diameter
    grid.N = n_sim;
    grid.N_ef = (packet_density/packet_mass)*pow(L,3)/grid.N;
    std::cout << "Each particle represents " << grid.N_ef << "packets of molecules/atoms\n";
    grid.num_dens = packet_density;
    grid.V = L*L*L;
    std::cout << "System volume (V): " << grid.V << std::endl;

    grid.cell_length = L / ncell;

    double v_init = sqrt(3*boltz*T/packet_mass);
    grid.v_mult = v_init;

    grid.a = 0.2;
    grid.v_max = 3*v_init; // max particle velocity
    grid.num_dt = 1e2;

    grid.create();
    grid.writeParticlesToDisk(out_dir + std::string("afterCreate-") + out_file);

    std::cout << "Particle width (d): " << grid.d << "\n";
    std::cout << "Mean free path: " << grid.lambda << std::endl; // should be 
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
