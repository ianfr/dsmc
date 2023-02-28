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
    std::string out_dir = "/Users/ian/Downloads/DSMC_OUT/";
    std::string cmd = "rm " + out_dir + "*.csv";
    std::system(cmd.c_str());
#endif
    std::string out_file = "particle.csv";

    int skip_every = 1; // skip every _ steps when writing out
    int n_zero = 10; // 10 digits for indexing output CSV series

    Grid grid;


    double boltz = 1.3806e-23; // J/K
    double mass = 6.63e-26; // mass argon
    double diam = 3.66e-10; // eff diam argon
    double T = 273; // temperature (K)
    double density = 1.78; // density of argon at STP (kg/m^3)
    double L = 100 * 1e-6; // x microns

    grid.d = diam; // kinetic particle diameter
    grid.N = 1e5;
    grid.N_ef = (density/mass)*pow(L,3)/grid.N;
    std::cout << "Each particle represents " << grid.N_ef << " molecules/atoms\n";
    grid.num_dens = density;
    grid.V = L*L*L;
    std::cout << "System volume (V): " << grid.V << std::endl;

    double ncell = 2; // number of cells in ONE dimension
    grid.cell_length = L / ncell;

    double v_init = sqrt(3*boltz*T/mass);
    grid.v_mult = v_init;

    grid.a = 0.2;
    grid.v_max = 3*v_init; // max particle velocity
    grid.num_dt = 1e2;

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
