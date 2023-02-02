//
// Created by Ian Friedrichs on 2/1/23.
//

#ifndef DSMC_CELL_H
#define DSMC_CELL_H

#include "Particle.h"

#include <vector>
#include <memory>
#include <random>
#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>

class Cell {
public:
    // Static variables, 'inline static' is C++17
    inline static float delta_t; // timestep
    inline static float v_max; // maximum particle velocity
    inline static int f_n; // number of molecules
    inline static float d; // particle diameter, uniform
    inline static float V_c; // cell volume, uniform for now
    inline static float v_mult; // initial velocity multiplier
    inline static float cell_length;

    // Data
    std::vector<std::shared_ptr<Particle>> m_part; // which members of the master list we're responsible for this iteration
    std::array<float,2> x_b; // x bounds: min,max
    std::array<float,2> y_b;
    std::array<float,2> z_b;

    // Methods
    void calculateCollisionsRejectionSampling(); // calculate new velocities but don't update
    void updatePositions(); // update the positions with euler (ok b/c no gravity for now)
    bool checkIfParticleInside(Particle &p); // check if a particle belongs in this box
    void initRandParticles(std::vector<std::shared_ptr<Particle>> &parts); // create random particles in this cell
    std::string str(); // string representation
};


#endif //DSMC_CELL_H
