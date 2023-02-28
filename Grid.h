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
    Cell* m_grid; // 3D Cell grid; m_grid(i,j,k) = m_grid[i * x * z + j * z + k]
    std::vector<Particle> m_part; // global list of all particles
    Vector3i grid_dims; // the number of cells in the x, y, and z directions
    int num_dt;
    int dim;
    double lambda; // mean free path
    double charlen; // characteristic length
    double mean_v; // mean velocity
    double a; // multiplier for timestep calculation
    double a_len; // multiplier for cell size
    int N; // the number of particles that are actually simulated
    double V; // system volume
    double num_dens; // number density; V = (N*N_ef)/num_dens; num_dens = N_phys_total / V


    // Cell static variables
    double delta_t; // timestep
    double v_max; // maximum particle velocity
    int N_ef; // the number of physical particles that each simulated particle represents
    double d; // particle diameter, uniform
    double cell_length; // length of one side of a (CUBIC) cell
    double v_mult; // initial velocity multiplier

    // Methods
    void calculateSystemVolume();
    void create();
    void calculateCollisionsRejectionSampling(); // calculate new velocities but don't update
    void updatePositions(); // update the positions with euler (ok b/c no gravity for now)
    void reassignParticlesToCells(); // brute force ensure all the particles are assigned to appropriate cells
    void enforceDomain(); // make sure particles don't leave the domain by reflecting them back in
    void writeParticlesToDisk(std::string filename);
    std::string str(); // get string representation of grid
    bool anyNullParticlePointers();

private:
    void setCellBoundaries();
    void addParticlesToCells();
};


#endif //DSMC_GRID_H
