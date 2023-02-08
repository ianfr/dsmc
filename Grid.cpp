//
// Created by Ian Friedrichs on 2/1/23.
//

#include "Grid.h"

void Grid::create() {
    // Allocate memory for the 3D grid
    m_grid = new Cell[grid_dims(0) * grid_dims(1) *grid_dims(2)];
    setCellBoundaries();

    // Allocate memory for the particles
    for (int i=0; i < f_n; i++) {
        Particle tmp;
        m_part.push_back(tmp);
    }

    // Set Cell static variables
    Cell::delta_t = delta_t;
    Cell::v_max = v_max;
    Cell::f_n = f_n;
    Cell::d = d;
    Cell::V_c = cell_length * cell_length * cell_length;
    Cell::cell_length = cell_length;
    Cell::v_mult = v_mult;

    addParticlesToCells();
}

// assume the grid starts in the "bottom left" at (0,0,0)
// WARNING: allocate grid first
void Grid::setCellBoundaries() {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);

    float c_x = 0;
    for (int i=0; i < x; i++) {
        float c_y = 0;
        for (int j=0; j < y; j++) {
            float c_z = 0;
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

    if (f_n % (x*y*z) != 0) {
        std::cout << "ERROR [Grid::addParticlesToCells] : the number of cells must evenly divide the number of particles " << std::endl;
        exit(1);
    }

    int num_p_per_cell = f_n / (x*y*z);
    int plist_idx = 0;
//    for (int i=0; i < x; i++) {
//        for (int j=0; j < y; j++) {
//            for (int k = 0; k < z; k++) {
//                std::vector<std::shared_ptr<Particle>> tmp_part(num_p_per_cell);
//                for (int idx = 0; idx < tmp_part.size(); idx++) {
//                    tmp_part[idx] = std::make_shared<Particle>(m_part[plist_idx + idx]);
//                }
//                m_grid[i * x * z + j * z + k].initRandParticles(tmp_part);
//                plist_idx += num_p_per_cell;
//            }
//        }
//    }
    for (int i = 0; i < x*y*z; i++) {
//        m_grid[i].m_part.reserve(num_p_per_cell);
        for (int idx = 0; idx < num_p_per_cell; idx++) {
            m_grid[i].m_part.push_back(&m_part[plist_idx + idx]);
            Vector3f r = Vector3f::Random() * cell_length;
            auto x_b = m_grid[i].x_b;
            auto y_b = m_grid[i].y_b;
            auto z_b = m_grid[i].z_b;
            Vector3f a = { ((x_b[0] + x_b[1]) / 2.0),
                           ((y_b[0] + y_b[1]) / 2.0),
                           ((z_b[0] + z_b[1]) / 2.0)}; // center of cell
            m_grid[i].m_part[idx]->pos = a + r;
            m_grid[i].m_part[idx]->vel = Vector3f::Random() * v_mult;
        }
        plist_idx += num_p_per_cell;
//        for (auto p : m_grid[i].m_part) std::cout << p->pos << std::endl;
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

void Grid::enforceDomain() {
    float eps = 1e-8;
    float x_min = 0;
    float y_min = 0;
    float z_min = 0;
    float x_max = grid_dims(0) * cell_length;
    float y_max = grid_dims(1) * cell_length;
    float z_max = grid_dims(2) * cell_length;
    for (auto &p : m_part) {
        float p_x = p.pos(0);
        float p_y = p.pos(1);
        float p_z = p.pos(2);

        if (p_x >= x_max) {
            p.pos(0) = x_max - eps;
            p.vel(0) *= -1;
        }
        if (p_y >= y_max) {
            p.pos(1) = y_max - eps;
            p.vel(1) *= -1;
        }
        if (p_z >= z_max) {
            p.pos(2) = z_max - eps;
            p.vel(2) *= -1;
        }

        if (p_x <= x_min) {
            p.pos(0) = x_min + eps;
            p.vel(0) *= -1;
        }
        if (p_y <= y_min) {
            p.pos(1) = y_min + eps;
            p.vel(1) *= -1;
        }
        if (p_z <= z_min) {
            p.pos(2) = z_min + eps;
            p.vel(2) *= -1;
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
