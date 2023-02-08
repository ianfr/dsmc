//
// Created by Ian Friedrichs on 2/1/23.
//

#ifndef DSMC_GRID_H
#define DSMC_GRID_H

#include "Cell.h"
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>

class Grid {
public:
    // Data
    // m_grid[i][j[k] = *(m_grid
    Cell* m_grid; // 3D Cell grid;
    std::vector<Particle> m_part; // global list of all particles
//    int n; // number of particles
    Vector3i grid_dims; // the number of cells in the x, y, and z directions

    // Cell static variables
    float delta_t; // timestep
    float v_max; // maximum particle velocity
    int f_n; // number of molecules
    float d; // particle diameter, uniform
    float cell_length; // length of one side of a (CUBIC) cell
    float v_mult; // initial velocity multiplier

    // Methods
    void create();
    void calculateCollisionsRejectionSampling(); // calculate new velocities but don't update
    void updatePositions(); // update the positions with euler (ok b/c no gravity for now)
    void reassignParticlesToCells(); // brute force ensure all the particles are assigned to appropriate cells
    void enforceDomain(); // make sure particles don't leave the domain by reflecting them back in
    void writeParticlesToDisk(std::string filename);
    std::string str(); // get string representation of grid
    bool anyNullParticlePointers();
    bool allParticlesInsideCells();

private:
    void setCellBoundaries();
    void addParticlesToCells();
};


#endif //DSMC_GRID_H
