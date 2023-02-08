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
    grid.f_n = 8*2; // number of particles
    grid.grid_dims = {2, 2, 2};
    grid.delta_t = 0.01;
    grid.v_max = 10.0; // max particle velocity
    grid.d = 0.5; // particle diameter
    grid.cell_length = 1.0;
    grid.v_mult = 1.0; // initial velocity multiplier

    grid.create();
    grid.writeParticlesToDisk(out_dir + std::string("afterCreate-") + out_file);

    // loop...
    grid.calculateCollisionsRejectionSampling();
    grid.writeParticlesToDisk(out_dir + std::string("afterCollisions-") + out_file);
    grid.updatePositions();
    grid.writeParticlesToDisk(out_dir + std::string("afterUpdatePositions-") + out_file);
    std::cout << "ANY NULL: " << grid.anyNullParticlePointers() << std::endl;

    grid.enforceDomain();
    grid.reassignParticlesToCells();
    std::cout << "ANY NULL: " << grid.anyNullParticlePointers() << std::endl;
    std::cout << "\n" << grid.str() << "\n";
//    grid.enforceDomain();
    std::cout << "\n" << grid.str() << "\n";
    grid.writeParticlesToDisk(out_dir + out_file);


    return 0;
}
