/* 
 @file      Boundaries.cpp
 @brief     contains functions for the Boundaries class
 @author    Luke Eure
 @date      January 9 2016
*/

#include "Boundaries.h"

/*
 @brief     constructor for Boundaries class
*/
Boundaries::Boundaries() {}

/*
 @brief     deconstructor
*/
Boundaries::~Boundaries() {}

/*
 @brief     set the surfaces into the geometry
 @paraam    axis 0, 1, or 2 corresponding to x y and z
 @param     side 0 or 1 corresponding to the minimum or maximum of the geometry
*/
void Boundaries::setSurface(Axes axis, min_max side, MCSurface* surface) {
    _surfaces[2*axis + side] = surface;
}

/*
 @brief     return the position of the surface
 @paraam    axis 0, 1, or 2 corresponding to x y and z
 @param     side 0 or 1 corresponding to the minimum or maximum of the geometry
 @return    the position of the surface within the geometry
*/
float Boundaries::getSurfaceCoord(int axis, int side) {
    return _surfaces[axis*2+side]->getPosition();
}

/*
 @brief     return the type of a surface
 @paraam    axis 0, 1, or 2 corresponding to x y and z
 @param     side 0 or 1 corresponding to the minimum or maximum of the geometry
 @return    the type of the surface, 0 = vacuum 1 = reflective
*/
BoundaryType Boundaries::getSurfaceType(int axis, int side) {
    return _surfaces[axis*2+side]->getType();
}

/*
 @brief     function that samples a random location within a bounding box.
 @details   a position along each axis is randomly and uniformally sampled in
            the bounding box. Each positin is set as the neutron's position
            along that axis.
 @param     neutron, the Neutron whose position will be set
*/
void Boundaries::sampleLocation(Neutron* neutron) {
    Point sampled_location;
    double width = getSurfaceCoord(0, MAX) - getSurfaceCoord(0, MIN);
    double coord = getSurfaceCoord(0, MIN) + width * neutron->arand();
    neutron->setPosition(0, coord);
    
    width = getSurfaceCoord(1, MAX) - getSurfaceCoord(1, MIN);
    coord = getSurfaceCoord(1, MIN) + width * neutron->arand();
    neutron->setPosition(1, coord);
    
    width = getSurfaceCoord(2, MAX) - getSurfaceCoord(2, MIN);
    coord = getSurfaceCoord(2, MIN) + width * neutron->arand();
    neutron->setPosition(2, coord);
}
