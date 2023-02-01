#include <iostream>

#include "Grid.h"

int main() {

    Grid grid;
    grid.f_n = 8*5; // number of particles
    grid.grid_dims = {2, 2, 2};
    grid.delta_t = 0.01;
    grid.v_max = 10.0; // max particle velocity
    grid.d = 0.5; // particle diameter
    grid.cell_length = 1.0;
    grid.v_mult = 1.0; // initial velocity multiplier

    grid.create();

    // loop...
    grid.calculateCollisionsRejectionSampling();
    grid.updatePositions();
    grid.reassignParticlesToCells();
    grid.enforceDomain();

    return 0;
}
