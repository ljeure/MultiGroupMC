/*
 @file      Mesh.h
 @brief     contains Mesh class
 @author    Luke Eure
 @date      January 9 2016
*/

#ifndef MESH_H
#define MESH_H

//#include <iostream>
#include <vector>
//#include <math.h>

#include "Boundaries.h"

class Mesh {
public:
    Mesh(Boundaries bounds, double delta_x, double delta_y, double delta_z,
            Material* default_material, int num_groups);
    virtual ~Mesh();

    void fluxAdd(std::vector <int> &cell, double distance, int group);
    void fluxClear();
    std::vector <std::vector <std::vector <std::vector <double> > > > getFlux();

private:

    /** the width of the cell along each axis */
    std::vector <double> _delta_axes;

    /** the minimum locations on the geometry in each direction */
    std::vector <double> _boundary_mins;

    /** the number of cells along each axis */
    std::vector <int> _axis_sizes;

    /** the neutron flux through each cell */
    std::vector <std::vector <std::vector <std::vector <double> > > > _flux;
    
    /** the number of energy groups */
    int _num_groups;
};

#endif
