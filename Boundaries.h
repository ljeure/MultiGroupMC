/* 
 @file      Boundaries.h
 @brief     container for Boundary objects
 @author    Luke Eure
 @date      January 6 2016
*/

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "Enumerations.h"
#include "Neutron.h"
#include "../../OpenMOC/src/Point.h"
#include "../../OpenMOC/src/Surface.h"

#include <vector>
#include <iostream>

class Boundaries {

public:
    Boundaries();
    virtual ~Boundaries();

    float getSurfaceCoord(int axis, int side);
    boundaryType getSurfaceType(int axis, int side);
    void setSurface(Axes axis, min_max side, double coord, boundaryType type);
    void sampleLocation(Neutron* neutron);

private:

    /** the coordinate of each geometry boundary along its axis*/
    double _surface_coords[6];

    /** the surface type of each geometry boundary */
    boundaryType _surface_types[6];
};

#endif
