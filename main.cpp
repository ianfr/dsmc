#include <iostream>

#include "Grid.h"

int main() {

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::string out_dir = R"(C:\Users\windowsuser\iCloudDrive\Projects\coding\monte-carlo\dsmc\OUT\)";
#else
    std::string out_dir = "/Users/ian/Library/Mobile Documents/com~apple~CloudDocs/Projects/coding/monte-carlo/dsmc/OUT/";
#endif
    std::string out_file = "particle.csv";

    int skip_every = 100; // skip every _ steps when writing out
    int n_zero = 10; // 10 digits for indexing

    Grid grid;
    grid.f_n = 8*5000; // number of particles
    grid.grid_dims = {2, 2, 2};
    grid.delta_t = 0.1;
    grid.num_dt = 10000;
    grid.v_max = 10.0; // max particle velocity
    grid.d = 0.001; // particle diameter
    grid.cell_length = 10.0;
    grid.v_mult = 1.0; // initial velocity multiplier

    grid.create();
    grid.writeParticlesToDisk(out_dir + std::string("afterCreate-") + out_file);

    int write_iter = 0;
    for (int iter=0; iter < grid.num_dt; iter++) {
        std::cout << "iter: " << iter << std::endl;
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
