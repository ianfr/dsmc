//
// Created by Ian Friedrichs on 2/1/23.
//

#include "Grid.h"

void Grid::create() {

    // Allocate memory for the particles
    for (int i=0; i < N; i++) {
        Particle tmp;
        m_part.push_back(tmp);
    }

    lambda = V / (sqrt(2) * N_ef * PI * d * d);
    charlen = std::pow((cell_length * V)/N_ef, 1.0/3.0);

    // Allocate memory for the 3D grid
    dim = (int) (pow(V, 1.0/3.0) / cell_length); // (domain side length / cell length)
    grid_dims = {dim, dim, dim};
    m_grid = new Cell[grid_dims(0) * grid_dims(1) *grid_dims(2)];

    // Set Cell static variables
    Cell::v_max = v_max;
    Cell::f_n = N_ef;
    Cell::d = d;
    Cell::V_c = cell_length * cell_length * cell_length;
    Cell::cell_length = cell_length;
    Cell::v_mult = v_mult;

    setCellBoundaries();
    addParticlesToCells();

    // Set the timestep
    std::vector<double> norms(N);
    std::transform(m_part.begin(), m_part.end(), norms.begin(),
                   [](auto &p) {return (p.vel).norm();});
    mean_v = std::accumulate(norms.begin(), norms.end(), 0.0f) / (double)N;
    delta_t = (a * cell_length) / mean_v;

    Cell::delta_t = delta_t;
}

// assume the grid starts in the "bottom left" at (0,0,0)
// WARNING: allocate grid first
void Grid::setCellBoundaries() {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);

    double c_x = 0;
    for (int i=0; i < x; i++) {
        double c_y = 0;
        for (int j=0; j < y; j++) {
            double c_z = 0;
            for (int k=0; k < z; k++) {
                m_grid[i * x * z + j * z + k].x_b = {c_x, c_x + cell_length};
                m_grid[i * x * z + j * z + k].y_b = {c_y, c_y + cell_length};
                m_grid[i * x * z + j * z + k].z_b = {c_z, c_z + cell_length};
                c_z += cell_length;
            }
            c_y += cell_length;
        }
        c_x += cell_length;
    }
#if PRINT_VERBOSE
    std::cout << "Added cell boundaries." << std::endl;
#endif
}

void Grid::addParticlesToCells() {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);

    int diff = 0;
    int num_p_per_cell;
    if (N % (x*y*z) != 0) {
        std::cout << "WARNING [Grid::addParticlesToCells] : the number of cells doesn't evenly divide the number of particles " << std::endl;
        diff = N % (x*y*z);
        num_p_per_cell = N / (x*y*z - 1);
    } else {
        num_p_per_cell = N / (x*y*z);
    }

    int plist_idx = 0;

    if (diff != 0) { // if it didn't divide perfectly
        // populate everything but the last cell (used for dumping the extra particles that didn't divide in properly)
        for (int i = 0; i < (x * y * z - 1); i++) {
            for (int idx = 0; idx < num_p_per_cell; idx++) {
                m_grid[i].m_part.push_back(&m_part[plist_idx + idx]);
                Vector3d r = Vector3d::Random() * cell_length;
                auto x_b = m_grid[i].x_b;
                auto y_b = m_grid[i].y_b;
                auto z_b = m_grid[i].z_b;
                Vector3d a = {((x_b[0] + x_b[1]) / 2.0),
                              ((y_b[0] + y_b[1]) / 2.0),
                              ((z_b[0] + z_b[1]) / 2.0)}; // center of cell
                m_grid[i].m_part[idx]->pos = a + r;
                m_grid[i].m_part[idx]->vel = Vector3d::Random().normalized() * v_mult;
            }
            plist_idx += num_p_per_cell;
        }

        int i = x * y * z - 1; // the last cell
        for (int idx = 0; idx < diff; idx++) {
            m_grid[i].m_part.push_back(&m_part[plist_idx + idx]);
            Vector3d r = Vector3d::Random() * cell_length;
            auto x_b = m_grid[i].x_b;
            auto y_b = m_grid[i].y_b;
            auto z_b = m_grid[i].z_b;
            Vector3d a = {((x_b[0] + x_b[1]) / 2.0),
                          ((y_b[0] + y_b[1]) / 2.0),
                          ((z_b[0] + z_b[1]) / 2.0)}; // center of cell
            m_grid[i].m_part[idx]->pos = a + r;
            m_grid[i].m_part[idx]->vel = Vector3d::Random().normalized() * v_mult;
        }
        plist_idx += num_p_per_cell;
    } else { // if it divides perfectly
        for (int i = 0; i < (x * y * z); i++) {
            for (int idx = 0; idx < num_p_per_cell; idx++) {
                m_grid[i].m_part.push_back(&m_part[plist_idx + idx]);
                Vector3d r = Vector3d::Random() * cell_length;
                auto x_b = m_grid[i].x_b;
                auto y_b = m_grid[i].y_b;
                auto z_b = m_grid[i].z_b;
                Vector3d a = {((x_b[0] + x_b[1]) / 2.0),
                              ((y_b[0] + y_b[1]) / 2.0),
                              ((z_b[0] + z_b[1]) / 2.0)}; // center of cell
                m_grid[i].m_part[idx]->pos = a + r;
                m_grid[i].m_part[idx]->vel = Vector3d::Random().normalized() * v_mult;
            }
            plist_idx += num_p_per_cell;
        }
    }



#if PRINT_VERBOSE
    std::cout << "\nDone adding particles to cells." << std::endl;
#endif
}

void Grid::calculateCollisionsRejectionSampling() {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);
    for (int i=0; i < x; i++) {
        for (int j=0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                m_grid[i * x * z + j * z + k].calculateCollisionsRejectionSampling();
            }
        }
    }
#if PRINT_VERBOSE
    std::cout << "Done calculating collisions." << std::endl;
#endif
}

void Grid::updatePositions() {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);
    for (int i=0; i < x; i++) {
        for (int j=0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                m_grid[i * x * z + j * z + k].updatePositions();
            }
        }
    }
#if PRINT_VERBOSE
    std::cout << "Done updating positions." << std::endl;
#endif
}

// helper function to check for overlap in Particle* lists
// adapted from https://stackoverflow.com/questions/19483663/vector-intersection-in-c
std::vector<Particle*> intersection(std::vector<Particle*> v1,
                                      std::vector<Particle*> v2){
    std::vector<Particle*> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}

void Grid::reassignParticlesToCells() {

    // have each Cell MOVE pointers for Particles that aren't inside them anymore
    std::vector<Particle*> need_to_move;

    // access the grid as a flat list
    // move
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);
    for (int i=0; i < x*y*z; i++) {
        for (int idx=0; idx < m_grid[i].m_part.size(); idx++) {
            if (!m_grid[i].checkIfParticleInside(m_grid[i].m_part[idx]->pos)) {
                need_to_move.push_back(m_grid[i].m_part[idx]);
                m_grid[i].m_part[idx] = nullptr;
            }
        }
    }

    // THEN use the erase-remove idiom to get rid of null pointers
    // https://stackoverflow.com/questions/11460810/finding-null-pointers-in-std-vectors
    // https://cplusplus.com/reference/memory/shared_ptr/operator%20bool/
    for (int i=0; i < x*y*z; i++) {
        m_grid[i].m_part.erase(std::remove_if(m_grid[i].m_part.begin(), m_grid[i].m_part.end(),
                                              [](auto ptr){
                                                        if (ptr == nullptr) {
                                                            return true;
                                                        }
                                                        return false;
                                                    }), m_grid[i].m_part.end());
    }

    for (int idx=0; idx < need_to_move.size(); idx++) {
        for (int i=0; i < x*y*z; i++) {
            if (need_to_move[idx] != nullptr) {
                if (m_grid[i].checkIfParticleInside(need_to_move[idx]->pos)) {
                    m_grid[i].m_part.push_back(need_to_move[idx]);
                    need_to_move[idx] = nullptr;
                }
            }
        }
    }

    // make sure there are no un-assigned particles
    need_to_move.erase(std::remove_if(need_to_move.begin(), need_to_move.end(),
                                [](auto ptr){
                                    if (ptr == nullptr) {
                                        return true;
                                    }
                                    return false;
                                }), need_to_move.end());
    if (need_to_move.size() > 0) {
        std::cout << "[ERROR]: unassigned particles: " << need_to_move.size() << "\n";
        for (auto p : need_to_move) {
            std::cout << p->pos << "\n";
        }
        exit(1);
    }

#if PRINT_VERBOSE
    std::cout << "Done reassigning particles to the appropriate cells." << std::endl;
#endif
}

// get a random unit vector in the hemisphere pointing away
    // from the minimum x
Vector3d getRandUnitHemVecMinX() {

    Vector3d a = Vector3d::Random().normalized();

    if (a.x() < 0) {
        a(0) *= -1;
    }
    return a;
}
// get a random unit vector in the hemisphere pointing away
// from the minimum y
Vector3d getRandUnitHemVecMinY() {

    Vector3d a = Vector3d::Random().normalized();

    if (a.y() < 0) {
        a(1) *= -1;
    }
    return a;
}
// get a random unit vector in the hemisphere pointing away
// from the minimum z
Vector3d getRandUnitHemVecMinZ() {

    Vector3d a = Vector3d::Random().normalized();

    if (a.z() < 0) {
        a(2) *= -1;
    }
    return a;
}
// get a random unit vector in the hemisphere pointing away
// from the maximum x
Vector3d getRandUnitHemVecMaxX() {

    Vector3d a = Vector3d::Random().normalized();

    if (a.x() > 0) {
        a(0) *= -1;
    }
    return a;
}
// get a random unit vector in the hemisphere pointing away
// from the maximum y
Vector3d getRandUnitHemVecMaxY() {

    Vector3d a = Vector3d::Random().normalized();

    if (a.y() > 0) {
        a(1) *= -1;
    }
    return a;
}
// get a random unit vector in the hemisphere pointing away
// from the maximum z
Vector3d getRandUnitHemVecMaxZ() {

    Vector3d a = Vector3d::Random().normalized();

    if (a.z() > 0) {
        a(2) *= -1;
    }
    return a;
}

void Grid::enforceDomain() {
    double eps = 1e-8;
    double x_min = 0;
    double y_min = 0;
    double z_min = 0;
    double x_max = grid_dims(0) * cell_length;
    double y_max = grid_dims(1) * cell_length;
    double z_max = grid_dims(2) * cell_length;
    for (auto &p : m_part) {
        double p_x = p.pos(0);
        double p_y = p.pos(1);
        double p_z = p.pos(2);

        if (p_x >= x_max) {
            p.pos(0) = x_max - eps;
            p.vel = p.vel.norm() * getRandUnitHemVecMaxX(); // diffuse reflection
        }
        if (p_y >= y_max) {
            p.pos(1) = y_max - eps;
            p.vel = p.vel.norm() * getRandUnitHemVecMaxY(); // diffuse reflection
        }
        if (p_z >= z_max) {
            p.pos(2) = z_max - eps;
            p.vel = p.vel.norm() * getRandUnitHemVecMaxZ(); // diffuse reflection
        }

        if (p_x <= x_min) {
            p.pos(0) = x_min + eps;
            p.vel = p.vel.norm() * getRandUnitHemVecMinX(); // diffuse reflection
//            p.vel(0) *= -0.9; // lose 10% energy when reflecting
        }

        if (p_y <= y_min) {
            p.pos(1) = y_min + eps;
            p.vel = p.vel.norm() * getRandUnitHemVecMinY(); // diffuse reflection
            p.vel(1) *= -1.05; // gain 5% energy when reflecting
        }
        if (p_z <= z_min) {
            p.pos(2) = z_min + eps;
            p.vel = p.vel.norm() * getRandUnitHemVecMinZ(); // diffuse reflection
        }

    }
#if PRINT_VERBOSE
    std::cout << "Done enforcing domain" << std::endl;
#endif
}

void Grid::writeParticlesToDisk(std::string filename) {
    std::ofstream outFile(filename, std::ios::trunc);
    outFile << "x, y, z, velmag\n";
    for (auto &p : m_part) {
        outFile << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << ", " << p.vel.norm() << "\n";
    }
}

std::string Grid::str() {
    std::ostringstream oss;
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);
    for (int i=0; i < x; i++) {
        for (int j=0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                oss << "CELL: (" << i << " " << j << " " << k << ")\n";
                oss << m_grid[i * x * z + j * z + k].str();
            }
        }
    }

    return oss.str();
}

bool Grid::anyNullParticlePointers() {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);
    for (int i=0; i < x*y*z; i++) {
        for (int idx=0; idx < m_grid[i].m_part.size(); idx++) {
            if ( ! m_grid[i].m_part[idx]) {
                return true;
            }
        }
    }
    return false;
}

void Grid::calculateSystemVolume() {
    V = N * (N_ef / num_dens);
}
