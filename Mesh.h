/*
 @file      Mesh.h
 @brief     contains Mesh class
 @author    Luke Eure
 @date      January 9 2016
*/

#ifndef MESH_H
#define MESH_H

#include <vector>

#include "Enumerations.h"
#include "Boundaries.h"

class Mesh {
public:
    Mesh(Boundaries bounds, int size_x, int size_y, int size_z,
            int num_groups);
    virtual ~Mesh();

    void fluxAdd(std::vector <int> &cell, double distance, int group);
    void fluxClear();
    std::vector <std::vector <std::vector <std::vector <double> > > > getFlux();

private:

    /** the number of cells along each axis */
    std::vector <int> _axis_sizes;

    /** the neutron flux through each cell */
    std::vector <std::vector <std::vector <std::vector <double> > > > _flux;
    
    /** the number of energy groups */
    int _num_groups;
};

#endif
