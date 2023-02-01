//
// Created by Ian Friedrichs on 2/1/23.
//

#include "Cell.h"

void Cell::calculateCollisionsRejectionSampling() {
    // number of candidates
    int n_c = (int)m_part.size();

    // calculate the number of candidates
    std::vector<float> norms(n_c);
    std::transform(m_part.begin(), m_part.end(), norms.begin(),
                   [](auto &p) {return (p->vel).norm();});
    float mean_v = std::accumulate(norms.begin(), norms.end(), 0.0f) / (float)n_c;
    int m_cand = (int) ((n_c * (n_c - 1) * f_n * M_PI * d*d * v_max * delta_t) / (2.0 * V_c));

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, (int) m_part.size() - 1); // define the range

    std::random_device rd01; // obtain a random number from hardware
    std::mt19937 gen01(rd01()); // seed the generator
    std::uniform_real_distribution<float> distr01(0, 1); // define the range

    for (int i=0; i < m_cand; i++) {
        // pick 2 particle indices
        int p1 = distr(gen);
        int p2 = distr(gen);

        // pick random proportion of maximum speed
        float v_tmp = distr01(gen01) * v_max;

        double v_r = (m_part[p1]->vel - m_part[p2]->vel).norm(); // relative speed
        if (v_r > v_tmp) {
            // calculate new velocities and update
            // http://corysimon.github.io/articles/uniformdistn-on-sphere/
            float theta = 2.0 * M_PI * distr01(gen01);
            float phi = acos(1 - 2 * distr01(gen01));
            // random unit vector
            Vector3f v_star = {
                    sin(theta) * cos(phi),
                    sin(theta) * sin(phi),
                    cos(theta)
            };

            // center of mass velocity
            Vector3f v_cm = 0.5 * (m_part[p1]->vel + m_part[p2]->vel);

            // update velocities
            m_part[p1]->vel = v_cm + 0.5 * v_star;
            m_part[p2]->vel = v_cm - 0.5 * v_star;
        }
    }



}

void Cell::updatePositions() {
    for (auto &p : m_part) {
        p->pos = p->pos + p->vel * delta_t;
    }
}

bool Cell::checkIfParticleInside(Particle &p) {

    if (p.pos(0) > x_b[0] && p.pos(0) < x_b[1]) {
        if (p.pos(1) > y_b[0] && p.pos(1) < y_b[1]) {
            if (p.pos(2) > z_b[0] && p.pos(2) < z_b[1]) {
                return true;
            }
        }
    }
    return false;
}

void Cell::initRandParticles(std::vector<std::shared_ptr<Particle>> &parts) {
    m_part.reserve(parts.size());
    m_part = parts;
    for (auto &p : parts) {
        // Random() generates random vectors with elements in [-1,1]
        Vector3f r = Vector3f::Random() * cell_length;
        Vector3f a = { ((x_b[0] + x_b[1]) / 2.0),
                       ((y_b[0] + y_b[1]) / 2.0),
                       ((z_b[0] + z_b[1]) / 2.0)}; // center of cell
        p->pos = a + r;
        p->vel = Vector3f::Random() * v_mult;
    }
}
