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
    for (int i=0; i < x; i++) {
        for (int j=0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                std::vector<std::shared_ptr<Particle>> tmp_part(num_p_per_cell);
                for (int idx = 0; idx < tmp_part.size(); idx++) {
                    tmp_part[idx] = std::make_shared<Particle>(m_part[plist_idx + idx]);
                }
                m_grid[i * x * z + j * z + k].initRandParticles(tmp_part);
                plist_idx += num_p_per_cell;
            }
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

// AWFUL way to do this, just for prototyping
void Grid::reassignParticlesToCells() {

    // Delete all the smart pointers that exist in Cells
//    int x = grid_dims(0);
//    int y = grid_dims(1);
//    int z = grid_dims(2);
//    for (int i=0; i < x; i++) {
//        for (int j=0; j < y; j++) {
//            for (int k = 0; k < z; k++) {
//                m_grid[i * x * z + j * z + k].m_part.clear();
//            }
//        }
//    }
//
//    for (int p_idx = 0; p_idx < m_part.size(); p_idx++) {
//        float p_x = m_part[p_idx].pos(0);
//        float p_y = m_part[p_idx].pos(1);
//        float p_z = m_part[p_idx].pos(2);
//
//        int x = grid_dims(0);
//        int y = grid_dims(1);
//        int z = grid_dims(2);
//        for (int i=0; i < x; i++) {
//            for (int j=0; j < y; j++) {
//                for (int k = 0; k < z; k++) {
//                    std::vector<std::shared_ptr<Particle>> tmp_part;
//                    if (m_grid[i * x * z + j * z + k].checkIfParticleInside(m_part[p_idx])) {
//                        auto tmp = std::make_shared<Particle>(m_part[p_idx]);
//                        m_grid[i * x * z + j * z + k].m_part.push_back(tmp);
//                    }
//                }
//            }
//        }
//    }

    // have each Cell MOVE smart pointers for Particles that aren't inside them anymore
    std::vector<std::shared_ptr<Particle>> need_to_move;

    // access the grid as a flat list
    // move
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);
    for (int i=0; i < x*y*z; i++) {
        for (int idx=0; idx < m_grid[i].m_part.size(); idx++) {
            if (!m_grid[i].checkIfParticleInside(*m_grid[i].m_part[idx].get())) {
                need_to_move.emplace_back(std::move(m_grid[i].m_part[idx]));
            }
        }
    }

    // THEN use the erase-remove idiom to get rid of null shared pointers
    // https://stackoverflow.com/questions/11460810/finding-null-pointers-in-std-vectors
    // https://cplusplus.com/reference/memory/shared_ptr/operator%20bool/
    for (int i=0; i < x*y*z; i++) {
        m_grid[i].m_part.erase(std::remove_if(m_grid[i].m_part.begin(), m_grid[i].m_part.end(),
                                              [](auto &ptr){ return (bool)ptr;}
                              ));
    }

//    int x = grid_dims(0);
//    int y = grid_dims(1);
//    int z = grid_dims(2);
//    for (int i=0; i < x; i++) {
//        for (int j=0; j < y; j++) {
//            for (int k = 0; k < z; k++) {
//                for (int idx=0; idx < m_grid[i * x * z + j * z + k].m_part.size(); idx++) {
//                    if (!m_grid[i * x * z + j * z + k].checkIfParticleInside(*m_grid[i * x * z + j * z + k].m_part[idx].get())) {
//                        need_to_move.emplace_back(std::move(m_grid[i * x * z + j * z + k].m_part[idx]));
//                        m_grid[i * x * z + j * z + k].m_part.erase(m_grid[i * x * z + j * z + k].m_part.begin() + idx);
//                    }
//                }
//            }
//        }
//    }

    // figure out where each of those smart pointers belongs and MOVE them to cells
    for (int i=0; i < x; i++) {
        for (int j=0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                for (int idx=0; idx < need_to_move.size(); idx++) {

                }
            }
        }
    }

#if PRINT_VERBOSE
    std::cout << "Done reassigning particles to the appropriate cells." << std::endl;
#endif
}

void Grid::enforceDomain() {
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

        if (p_x > x_max) {
            p.pos(0) = x_max;
            p.vel(0) *= -1;
        }
        if (p_y > y_max) {
            p.pos(1) = y_max;
            p.vel(1) *= -1;
        }
        if (p_z > z_max) {
            p.pos(2) = z_max;
            p.vel(2) *= -1;
        }

        if (p_x < x_min) {
            p.pos(0) = x_min;
            p.vel(0) *= -1;
        }
        if (p_y < y_min) {
            p.pos(1) = y_min;
            p.vel(1) *= -1;
        }
        if (p_z < z_min) {
            p.pos(2) = z_min;
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
