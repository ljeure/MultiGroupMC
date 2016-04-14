/*
 @file      Mesh.cpp
 @brief     contains functions for the Mesh class
 @author    Luke Eure
 @date      January 9 2016
*/

#include "Mesh.h"

/*
 @brief     constructor for Mesh class
*/
Mesh::Mesh(Boundaries bounds, double delta_x, double delta_y, double delta_z,
        Material* default_material, int num_groups) {

    // save deltas 
    _delta_axes.push_back(delta_x);
    _delta_axes.push_back(delta_y);
    _delta_axes.push_back(delta_z);

    // save number of groups
    _num_groups = num_groups;

    // save boundary mins
    for (int axis=0; axis<3; ++axis) {
        _boundary_mins.push_back(bounds.getSurfaceCoord(axis, 0));
    }

    // save axis sizes
    for (int axis=0; axis<3; ++axis) {
        int size = (bounds.getSurfaceCoord(axis, MAX) - _boundary_mins[axis])
            / _delta_axes[axis];
        _axis_sizes.push_back(size);
    }
    
    // resize _flux and set all its elements = 0
    _flux.resize(_num_groups);
    for (int g=0; g<_num_groups; ++g) {
        _flux[g].resize(_axis_sizes[0]);
        for (int i=0; i<_axis_sizes[0]; ++i) {
            _flux[g][i].resize(_axis_sizes[1]);
            for (int j=0; j<_axis_sizes[1]; ++j) {
                _flux[g][i][j].resize(_axis_sizes[2]);
                for (int k=0; k<_axis_sizes[2]; ++k) {
                    _flux[g][i][j][k] = 0.0;
                }
            }
        }
    }
}

/*
 @brief     deconstructor for Mesh class
*/
Mesh::~Mesh() {}


/*
 @brief     add the distance a neutron has traveled within the cell to the flux
            array
 @param     cell a vector containing a cell
 @param     distance a distance to be added to the cell flux
 @param     group a group to which this distance should be added
*/
void Mesh::fluxAdd(std::vector <int> &cell, double distance, int group) {
    _flux[group][cell[0]][cell[1]][cell[2]] += distance;
}

/*
 @brief     set the value of each element in the flux array to 0
*/
void Mesh::fluxClear() {
    for (int g=0; g<_num_groups; ++g) {
        for (int i=0; i<_axis_sizes[0]; ++i) {
            for (int j=0; j<_axis_sizes[1]; ++j) {
                for (int k=0; k<_axis_sizes[2]; ++k) {
                    _flux[g][i][j][k] = 0.0;
                }
            }
        }
    }
}

/*
 @brief     return the flux array
 @return    returns the 4d flux vector
*/
std::vector <std::vector <std::vector <std::vector <double> > > > 
        Mesh::getFlux() {
    return _flux;
}

/*
   searching among files grep -r "getFSRId" OpenMOC

   in montecarlo
    initlialize fsrs function
    loop that creates points in each region
    findFSRId(point) for each point. This adds the FSR to the geometry
   


   index cells by fsr_id
   cell.getFSRId(cell)

   each surface splits into half spaces
    
   */
