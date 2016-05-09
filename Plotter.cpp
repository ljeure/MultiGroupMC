/*
 @file      Plotter.cpp
 @brief     plotting function
 @author    Luke Eure
 @date      January 28 2016
*/

#include "Plotter.h"

void printFluxToFile(std::vector <double> &flux_fsr) {
    std::ofstream out ("flux_fsr_plot.txt");

    for (int i=0; i<flux_fsr.size(); ++i) {
        out << flux_fsr[i] << " \n";
    }
}
