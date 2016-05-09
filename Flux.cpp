/*
 @file      Flux.cpp
 @brief     contains functions for the Flux class
 @author    Luke Eure
 @date      January 9 2016
*/

#include "Flux.h"

/*
 @brief     constructor for Flux class
*/
Flux::Flux(int num_fsr, int num_groups) {

    // save number of groups
    _num_groups = num_groups;

    // save number of elements in the flux array
    _num_fsr = num_fsr;

    // resize _fsr_flux and set its elements = 0
    _flux.resize(_num_groups*_num_fsr);
    for (int i=0; i<_num_groups*_num_fsr; ++i)
        _flux[i] = 0.0;

}

/*
 @brief     deconstructor for Flux class
*/
Flux::~Flux() {}

/*
 @brief     add the distance a neutron has traveled within the cell to the flux
            array
 @param     group the energy group for the distance to be added to
 @param     fsr the id of the flat source region for the distance to be added to
 @param     distance a distance to be added to the cell flux
*/
void Flux::add(int group, int fsr, double distance) {
    _flux[group*_num_fsr+fsr] += distance;
}

/*
 @brief     set the value of each element in the flux array to 0
*/
void Flux::clear() {
    for (int i=0; i<_num_groups*_num_fsr; ++i)
        _flux[i] = 0.0;
}

/*
 @brief     return the fsr flux array
 @return    returns the 1d flux vector
*/
std::vector <double> Flux::getFlux() {
    return _flux;
}
