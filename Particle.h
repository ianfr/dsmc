//
// Created by Ian Friedrichs on 2/1/23.
//

#ifndef DSMC_PARTICLE_H
#define DSMC_PARTICLE_H

#include "Eigen3/Eigen/Dense"

#define PRINT_VERBOSE 1

using namespace Eigen;

struct Particle {
    Vector3f pos = {0,0,0}; // position in R3
    Vector3f vel = {0,0,0}; // velocity in R3
};


#endif //DSMC_PARTICLE_H
