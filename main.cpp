/* 
 @file      Main.cpp
 @brief     creates geometry and materials to run a Monte Carlo simulation
 @author    Luke Eure
 @date      January 9 2016
*/

#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "Tally.h"
#include "Flux.h"
#include "Plotter.h"
#include "Neutron.h"
#include "Monte_carlo_solver.h"
#include "../../OpenMOC/src/Universe.h"
#include "../../OpenMOC/src/Material.h"
#include "../../OpenMOC/src/Cell.h"
#include "../../OpenMOC/src/Universe.h"
#include "../../OpenMOC/src/Geometry.h"

int main() {

    // create openmoc surfaces and set their boundary types
    XPlane* x_min = new XPlane(-2.0, 0, "x_min");
    XPlane* x_max = new XPlane(2.0, 1, "x_max");
    YPlane* y_min = new YPlane(-2.0, 2, "y_min");
    YPlane* y_max = new YPlane(2.0, 3, "y_max");
    x_min->setBoundaryType(VACUUM);
    x_max->setBoundaryType(VACUUM);
    y_min->setBoundaryType(VACUUM);
    y_max->setBoundaryType(VACUUM);

    // number of energy groups
    int num_groups = 2;
    
    // create fuel
    Material* fuel = new Material(1, "fuel");
    fuel->setNumEnergyGroups(num_groups);
    fuel->setSigmaTByGroup(2.0/9.0, 1);
    fuel->setSigmaTByGroup(5.0/6.0, 2);
    fuel->setSigmaFByGroup(1.0/480.0, 1);
    fuel->setSigmaFByGroup(1.0/16.0, 2);
    fuel->setNuSigmaFByGroup(2.4/480.0, 1);
    fuel->setNuSigmaFByGroup(2.4/16.0, 2);
    fuel->setSigmaSByGroup(71.0/360.0, 1, 1);
    fuel->setSigmaSByGroup(.02, 1, 2);
    fuel->setSigmaSByGroup(0.0, 2, 1);
    fuel->setSigmaSByGroup(11.0/15.0, 2, 2);
    fuel->setChiByGroup(1.0, 1);
    fuel->setChiByGroup(0.0, 2);

    // create moderator
    Material* moderator = new Material(2, "moderator");
    moderator->setNumEnergyGroups(num_groups);
    moderator->setSigmaTByGroup(2.0/9.0, 1);
    moderator->setSigmaTByGroup(5.0/3.0, 2);
    moderator->setSigmaFByGroup(0.0, 1);
    moderator->setSigmaFByGroup(0.0, 2);
    moderator->setNuSigmaFByGroup(0.0, 1);
    moderator->setNuSigmaFByGroup(0.0, 2);
    moderator->setSigmaSByGroup(71.0/360.0, 1, 1);
    moderator->setSigmaSByGroup(.025, 1, 2);
    moderator->setSigmaSByGroup(0.0, 2, 1);
    moderator->setSigmaSByGroup(47.0/30.0, 2, 2);
    moderator->setChiByGroup(0.0, 1);
    moderator->setChiByGroup(0.0, 2);

    // create cells
    Cell* root_cell = new Cell(0, "root");
    root_cell->addSurface(1, x_min);
    root_cell->addSurface(-1, x_max);
    root_cell->addSurface(1, y_min);
    root_cell->addSurface(-1, y_max);
    
    Cell* moderator_cell = new Cell(1, "moderator");
    moderator_cell->setFill(moderator);
    
    Cell* fuel_cell = new Cell(2, "fuel");
    fuel_cell->setFill(fuel);

    // create universes
    Universe* root_universe = new Universe(0, "root universe");
    root_universe->addCell(root_cell);

    Universe* moderator_universe = new Universe(1, "moderator universe");
    moderator_universe->addCell(moderator_cell);

    Universe* fuel_universe = new Universe(2, "fuel universe");
    fuel_universe->addCell(fuel_cell);

    // create lattice
    int numXLat = 9;
    int numYLat = 9;
    int numZLat = 1;
    Lattice* lattice = new Lattice();
    lattice->setWidth(4.0/9.0, 4.0/9.0);

    // create universe array for input into lattice
    Universe* universe_array [numXLat*numYLat];
    for (int i=0; i<numXLat; ++i) {
        for (int j=0; j<numYLat; ++j) {
            universe_array[j*numXLat+i] = moderator_universe;
        }
    }
    for (int i=3; i<6; ++i) {
        for (int j=3; j<6; ++j) {
            universe_array[j*numXLat+i] = fuel_universe;
        }
    }
    Universe** universe_array_pointer;
    universe_array_pointer = universe_array;

    // set universes in lattice
    lattice->setUniverses(numZLat, numYLat, numXLat, universe_array_pointer);

    // fill root cell with lattice
    root_cell->setFill(lattice);

    // create geometry
    Geometry* geometry = new Geometry();
    geometry->setRootUniverse(root_universe);
    
//----------------------------------------------------------------------------//
   
    // simulate neutron histories
    int num_neutrons = 10000;
    int num_batches = 3;
    
    // create solver
    MCSolver solver;
    solver.setGeometry(geometry);
    solver.setRootCell(root_cell);

    // initialize solver function
    solver.initializeFSRs(lattice);
    solver.initializeFlux();

    // simulate neutrons
    solver.computeEigenValue(num_neutrons, num_batches, num_groups);
    std::cout << "gets here\n";

    // plot neutron flux
    std::vector <double> flux = solver.getFlux()->getFlux();
    printFluxToFile(flux);

    // run python script to get flux plots
    system("python Flux_parser.py");

    std::cout << std::endl;

    delete geometry;
    delete x_min;
    delete x_max;
    delete y_min;
    delete y_max;
/*  
    get seg faults when I try to delete these
    delete lattice;
    delete fuel_cell;
    delete moderator_cell;
    delete fuel_universe;
    delete moderator_universe;
    delete root_universe;
    delete moderator;
    delete fuel;
*/

    return 0;
}


/*

main should look like profile/models homogeneous

   */


