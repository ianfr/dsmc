#include <iostream>

#include "Grid.h"

int main() {

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::string out_dir = R"(C:\Users\windowsuser\iCloudDrive\Projects\coding\monte-carlo\dsmc\OUT\)";
#else
    std::string out_dir = "/Users/ian/Library/Mobile Documents/com~apple~CloudDocs/Projects/coding/monte-carlo/dsmc/OUT/";
#endif
    std::string out_file = "particle.csv";

    Grid grid;
    grid.f_n = 8*5000; // number of particles
    grid.grid_dims = {2, 2, 2};
    grid.delta_t = 0.1;
    grid.num_dt = 100;
    grid.v_max = 10.0; // max particle velocity
    grid.d = 0.001; // particle diameter
    grid.cell_length = 10.0;
    grid.v_mult = 1.0; // initial velocity multiplier

    grid.create();
    grid.writeParticlesToDisk(out_dir + std::string("afterCreate-") + out_file);

    for (int iter=0; iter < grid.num_dt; iter++) {
        std::cout << "iter: " << iter << std::endl;
        grid.calculateCollisionsRejectionSampling();
        //    grid.writeParticlesToDisk(out_dir + std::string("afterCollisions-") + out_file);
        grid.updatePositions();
        //    grid.writeParticlesToDisk(out_dir + std::string("afterUpdatePositions-") + out_file);
        grid.enforceDomain();
        grid.reassignParticlesToCells();
//        std::cout << "\n" << grid.str() << "\n";
        std::string iter_str = std::to_string(iter);
        grid.writeParticlesToDisk(out_dir + iter_str + std::string("-") + out_file);
    }

    return 0;
}
