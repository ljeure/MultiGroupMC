/*
 @file      Flux.h
 @brief     contains Flux class
 @author    Luke Eure
 @date      January 9 2016
*/

#ifndef Flux_H
#define Flux_H

#include <vector>

#include "Enumerations.h"

class Flux {
public:
    Flux(int num_fsr, int num_groups);
    virtual ~Flux();

    void add(int group, int fsr, double distance);
    void clear();
    std::vector <double> getFlux();

private:
    
    /** the neutron flux through each flat source region */
    std::vector <double> _flux;

    /** the number of energy groups */
    int _num_groups;

    /** the number of FSRs in the geometry */
    int _num_fsr;
};

#endif
