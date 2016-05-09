/* 
 @file      Enumerations.h
 @brief     enumerations used in Monte Carlo
 @author    Luke Eure
 @date      April 14 2016
*/

#ifndef ENUMERATIONS_H
#define ENUMERATIONS_H

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

#endif
