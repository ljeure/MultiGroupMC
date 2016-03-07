/* 
 @file      Boundary.h
 @brief     holds the Boundary class
 @author    Luke Eure
 @date      February 13 2016
*/

#ifndef Surface_H
#define Surface_H

//#include "precision.h"

enum BoundaryType {
    MCVACUUM,
    MCREFLECTIVE,
    MCBOUNDARY_NONE
};

// enumerations used when talking about boundaries
enum min_max {MIN, MAX};

enum Axes {
    X,
    Y,
    Z
};

class MCSurface {

public:
    MCSurface(BoundaryType type, double position);
    virtual ~MCSurface();

    double getPosition();
    BoundaryType getType();

private:

    /** the type of the boundary, vacuum, reflective, or none */
    BoundaryType _boundary_type;

    /** the position of the boundary */
    double _position;
};


#endif
