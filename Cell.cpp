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

bool Cell::checkIfParticleInside(Vector3f pos) {

    if (pos(0) > x_b[0] && pos(0) < x_b[1]) {
        if (pos(1) > y_b[0] && pos(1) < y_b[1]) {
            if (pos(2) > z_b[0] && pos(2) < z_b[1]) {
                return true;
            }
        }
    }
    return false;
}

// not modifying the particles for some reason
//void Cell::initRandParticles(std::vector<std::shared_ptr<Particle>> &parts) {
////void Cell::initRandParticles(std::vector<std::shared_ptr<Particle>> parts) {
//    m_part.reserve(parts.size());
////    m_part = parts;
//    m_part = std::move(parts);
////    for (auto &p : parts) {
//    for (auto &p : m_part) {
//        // Random() generates random vectors with elements in [-1,1]
//        Vector3f r = Vector3f::Random() * cell_length;
//        Vector3f a = { ((x_b[0] + x_b[1]) / 2.0),
//                       ((y_b[0] + y_b[1]) / 2.0),
//                       ((z_b[0] + z_b[1]) / 2.0)}; // center of cell
//        p->pos = a + r;
//        p->vel = Vector3f::Random() * v_mult;
//        std::cout << p->pos << std::endl;
//    }
//}

std::string Cell::str() {

    std::ostringstream oss;

    oss << "x_bounds: " << x_b[0] << " " << x_b[1] << "\n";
    oss << "y_bounds: " << y_b[0] << " " << y_b[1] << "\n";
    oss << "z_bounds: " << z_b[0] << " " << z_b[1] << "\n";

    for (auto p : m_part) {
        oss << "pos: " << p->pos(0) << " " << p->pos(1) << " " << p->pos(2)
           << " vel: " << p->vel(0) << " " << p->vel(1) << " " << p->vel(2) << "\n";
    }

    return oss.str();
}
