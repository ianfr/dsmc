//
// Created by Ian Friedrichs on 2/1/23.
//

#ifndef DSMC_GRID_H
#define DSMC_GRID_H

#include "Cell.h"

class Grid {
public:
    // Data
    // m_grid[i][j[k] = *(m_grid
    Cell* m_grid; // 3D Cell grid;
    std::vector<Particle> m_part; // global list of all particles
    int n; // number of particles
    Vector3i grid_dims; // the number of cells in the x, y, and z directions

    // Cell static variables
    float delta_t; // timestep
    float v_max; // maximum particle velocity
    int f_n; // number of molecules
    float d; // particle diameter, uniform
    float cell_length; // length of one side of a (CUBIC) cell
    float v_mult; // initial velocity multiplier

    void create();

private:
    void setCellBoundaries();
    std::shared_ptr<Cell> getCellPtr(int i, int j, int k);
    void addParticlesToCells();
};


#endif //DSMC_GRID_H
