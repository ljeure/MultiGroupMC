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

#include "Boundaries.h"
#include "Tally.h"
#include "Flux.h"
#include "Plotter.h"
#include "Neutron.h"
#include "Monte_carlo.h"
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
//    ZPlane* z_min = new ZPlane(-2.0, 4, "z_min");
//    ZPlane* z_max = new ZPlane(2.0, 5, "z_max");
    x_min->setBoundaryType(VACUUM);
    x_max->setBoundaryType(VACUUM);
    y_min->setBoundaryType(VACUUM);
    y_max->setBoundaryType(VACUUM);
//    z_min->setBoundaryType(REFLECTIVE);
//    z_max->setBoundaryType(REFLECTIVE);

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
//    root_cell->addSurface(1, z_min);
//    root_cell->addSurface(-1, z_max);
    
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
    
    // add FSR points to geometry
    for (int i=0; i<lattice->getNumX(); ++i) {
        double xPos = lattice->getMinX() + lattice->getWidthX()*(i + .5);
        
        for (int j=0; j<lattice->getNumX(); ++j) {
            double yPos = lattice->getMinY() + lattice->getWidthY()*(j + .5);
        
            for (int k=0; k<lattice->getNumZ(); ++k) {
                double zPos = lattice->getMinZ()
                    + lattice->getWidthZ()*(k + .5);
                LocalCoords* localCoordFSR = new LocalCoords(xPos, yPos, zPos);
                localCoordFSR->setUniverse(root_universe);
                geometry->findFirstCell(localCoordFSR);
                
                int fsr = geometry->findFSRId(localCoordFSR);

                delete localCoordFSR;
            }
        }
    }
    geometry->initializeFSRVectors();
    
//----------------------------------------------------------------------------//
    
    // create geometry with surfaces
    // z boundaries are automatically reflective because OpenMOC doesn't
    // support 3D
    Boundaries test_boundary;
    test_boundary.setSurface(X, MAX, root_cell->getMaxX(),
            root_cell->getMaxXBoundaryType());
    test_boundary.setSurface(X, MIN, root_cell->getMinX(),
            root_cell->getMinXBoundaryType());
    test_boundary.setSurface(Y, MAX, root_cell->getMaxY(),
            root_cell->getMaxYBoundaryType());
    test_boundary.setSurface(Y, MIN, root_cell->getMinY(),
            root_cell->getMinYBoundaryType());
//    test_boundary.setSurface(Z, MAX, root_cell->getMaxZ(), REFLECTIVE);
//    test_boundary.setSurface(Z, MIN, root_cell->getMinZ(), REFLECTIVE);

    // create flux
    Flux test_flux(geometry->getNumFSRs(), fuel->getNumEnergyGroups());

//----------------------------------------------------------------------------//
    
    // simulate neutron histories
    int num_neutrons = 1000;
    int num_batches = 3;
    
    generateNeutronHistories(num_neutrons, test_boundary, test_flux, lattice,
            num_batches, num_groups, geometry, root_universe);
    
    // plot neutron flux
    std::vector <double> flux= test_flux.getFlux();
    printFluxToFile(flux);

    // run python script to get flux plots
    system("python Flux_parser.py");

    std::cout << std::endl;

    delete geometry;
    delete x_min;
    delete x_max;
    delete y_min;
    delete y_max;
//    delete z_min;
//    delete z_max;
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
   make monte carlo a solver
   should have same inputs as cpu solver + num_neutrons

   cpu has comput.eval
   MC should have compute eval (num_inters)

   look at cpusolver.h/.cpp, solver.h/.cpp


    
    cpu solver is sublcass of solver
    monte carlo solver is also a sublcass of solver

   for solver:
   if it has memberFunction() = 0;
    I need memberFunction( same parameters) {
         return 0 or NULL or whatever;}

*/

