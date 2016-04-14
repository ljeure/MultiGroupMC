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
Mesh::Mesh(Boundaries bounds, int size_x, int size_y, int size_z,
        int num_groups) {
    
    // save axis sizes
    _axis_sizes.push_back(size_x);
    _axis_sizes.push_back(size_y);
    _axis_sizes.push_back(size_z);

    // save number of groups
    _num_groups = num_groups;

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
   index cells by fsr_id
   cell.getFSRId(cell)

   each surface splits into half spaces
    
   */
