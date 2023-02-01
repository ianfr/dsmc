//
// Created by Ian Friedrichs on 2/1/23.
//

#include "Grid.h"

void Grid::create() {
    // Allocate memory for the 3D grid
    m_grid = new Cell[grid_dims(0) * grid_dims(1) *grid_dims(2)];
    setCellBoundaries();

    // Allocate memory for the particles
    for (int i=0; i < n; i++) {
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
                auto cell = getCellPtr(i, j, k);
                cell->x_b = {c_x, c_x + cell_length};
                cell->y_b = {c_y, c_y + cell_length};
                cell->z_b = {c_z, c_z + cell_length};
                c_z += cell_length;
            }
            c_y += cell_length;
        }
        c_x += cell_length;
    }
}

// row-major, see https://stackoverflow.com/questions/10151418/indexing-a-3-dimensional-array-using-a-single-contiguous-block-of-memory
std::shared_ptr<Cell> Grid::getCellPtr(int i, int j, int k) {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);

    auto tmp = &m_grid[i * x * z + j * z + k];

    return std::shared_ptr<Cell>(tmp);
}

void Grid::addParticlesToCells() {
    int x = grid_dims(0);
    int y = grid_dims(1);
    int z = grid_dims(2);

    // CONTINUE HERE...

    for (int i=0; i < x; i++) {
        for (int j=0; j < y; j++) {
            for (int k=0; k < z; k++) {
                auto cell = getCellPtr(i, j, k);

            }
        }
    }
}
