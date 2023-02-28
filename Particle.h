//
// Created by Ian Friedrichs on 2/1/23.
//

#ifndef DSMC_PARTICLE_H
#define DSMC_PARTICLE_H

#define PI 3.14159265

#include "Eigen3/Eigen/Dense"

#define PRINT_VERBOSE 0

using namespace Eigen;

struct Particle {
    Vector3d pos = {0,0,0}; // position in R3
    Vector3d vel = {0,0,0}; // velocity in R3
};


#endif //DSMC_PARTICLE_H
