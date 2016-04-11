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

#include "Surface.h"
#include "Boundaries.h"
#include "Tally.h"
#include "Mesh.h"
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
    XPlane* x_max = new XPlane(2.1, 1, "x_max");
    YPlane* y_min = new YPlane(-2.0, 2, "y_min");
    YPlane* y_max = new YPlane(2.0, 3, "y_max");
    ZPlane* z_min = new ZPlane(-2.0, 4, "z_min");
    ZPlane* z_max = new ZPlane(2.0, 5, "z_max");
    x_min->setBoundaryType(VACUUM);
    x_max->setBoundaryType(VACUUM);
    y_min->setBoundaryType(VACUUM);
    y_max->setBoundaryType(VACUUM);
    z_min->setBoundaryType(REFLECTIVE);
    z_max->setBoundaryType(REFLECTIVE);

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
    Material* moderator = new Material(0, "moderator");
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
    root_cell->addSurface(1, z_min);
    root_cell->addSurface(-1, z_max);
    
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
    Lattice lattice;
    lattice.setWidth(4.0/9.0, 4.0/9.0);
    
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
    lattice.setUniverses(numZLat, numYLat, numXLat, universe_array_pointer);

    // create geometry
    Geometry* geometry = new Geometry();
    geometry->setRootUniverse(root_universe);

//-------------------------------------------------------------------------
    // create geometry with surfaces
    Boundaries test_boundary;
    test_boundary.setSurface(X, MAX, x_max);
    test_boundary.setSurface(X, MIN, x_min);
    test_boundary.setSurface(Y, MAX, y_max);
    test_boundary.setSurface(Y, MIN, y_min);
    test_boundary.setSurface(Z, MAX, z_max);
    test_boundary.setSurface(Z, MIN, z_min);

    // create array with z_coordinate max and mins for 2D OpenMOC
    std::vector <double> z_bounds (2);
    z_bounds[0] = test_boundary.getSurfaceCoord(2, 0);
    z_bounds[1] = test_boundary.getSurfaceCoord(2, 1);

    // create mesh
    Mesh test_mesh(test_boundary, 4.0/9.0, 4.0/9.0, 4.0, moderator,
            num_groups);

    // fill mesh with some material
    static const double a_fuel_limits [6] =
    { -2.0/3.0,     2.0/3.0,
      -2.0/3.0,     2.0/3.0,
      -2.0,         2.0 };
    std::vector <std::vector <double> > fuel_limits (3,
            (std::vector <double> (2)));
    for (int i=0; i<3; ++i) {
        for (int j=0; j<2; ++j) {
            fuel_limits[i][j] = a_fuel_limits[i*2+j];
        }
    }
    test_mesh.fillMaterials(fuel, fuel_limits);

    // test fsrid
    std::cout << "testing fsrid\n";
    LocalCoords* test = new LocalCoords(.1,.1,.1);
    test->setUniverse(root_universe);
    //test->setCell(root_cell);
    std::cout << "testing fsrid\n";
    std::cout << "fsr test " << geometry->findFSRId(test) << std::endl;

    // add FSR points to geometry
    for (int i=0; i<lattice.getNumX(); ++i) {
        double xPos = lattice.getMinX() + lattice.getWidthX()*(i + 1/2);
        
        for (int j=0; j<lattice.getNumX(); ++j) {
            double yPos = lattice.getMinY() + lattice.getWidthY()*(j + 1/2);
        
            for (int k=0; k<lattice.getNumZ(); ++k) {
                double zPos = lattice.getMinZ() + lattice.getWidthZ()*(k + 1/2);
                LocalCoords* rootCoordFSR = new LocalCoords(xPos, yPos, zPos);
                Universe* cellUniverse = lattice.getUniverse(i, j, k);
                Cell* cellForFSR = root_universe->findCell(
                        rootCoordFSR);
                LocalCoords* localCoordFSR = new LocalCoords(0, 0, 0);
                localCoordFSR->setUniverse(cellUniverse);
                localCoordFSR->setCell(cellForFSR);
                
                std::cout << "fsrid " << i << " " << j << " " << k
                    << ": " << geometry->findFSRId(localCoordFSR) << std::endl;

            }
        }
    }


    /*
    // simulate neutron histories
    int num_neutrons = 990;
    int num_batches = 3;
    
    generateNeutronHistories(num_neutrons, test_boundary,
            test_mesh, lattice, num_batches, num_groups, geometry, z_bounds);
    
    // plot neutron flux
    std::vector <std::vector <std::vector <std::vector <double> > > > flux =
       test_mesh.getFlux();
    printFluxToFile(flux);

    // run python script to get flux plots
    system("python Flux_parser.py");
*/

    std::cout << std::endl;
    return 0;

    }
