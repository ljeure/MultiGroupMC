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
//#include "../../OpenMOC/src/Surface.cpp"

int main() {

    // create openmoc surfaces and set their boundary types
    XPlane x_min(-2.0, 0, "x_min");
    XPlane x_max(2.1, 1, "x_max");
    YPlane y_min(-2.0, 2, "y_min");
    YPlane y_max(2.0, 3, "y_max");
    ZPlane z_min(-2.0, 4, "z_min");
    ZPlane z_max(2.0, 5, "z_max");
    x_min.setBoundaryType(VACUUM);
    x_max.setBoundaryType(VACUUM);
    y_min.setBoundaryType(VACUUM);
    y_max.setBoundaryType(VACUUM);
    z_min.setBoundaryType(REFLECTIVE);
    z_max.setBoundaryType(REFLECTIVE);

    // create surface pointers
    Surface* x_min_point = &x_min;
    Surface* x_max_point = &x_max;
    Surface* y_min_point = &y_min;
    Surface* y_max_point = &y_max;
    Surface* z_min_point = &z_min;
    Surface* z_max_point = &z_max;

    // number of energy groups
    int num_groups = 2;
    
    // create fuel
    Material fuel(1, "fuel");
    fuel.setNumEnergyGroups(num_groups);
    fuel.setSigmaTByGroup(2.0/9.0, 1);
    fuel.setSigmaTByGroup(5.0/6.0, 2);
    fuel.setSigmaFByGroup(1.0/480.0, 1);
    fuel.setSigmaFByGroup(1.0/16.0, 2);
    fuel.setNuSigmaFByGroup(2.4/480.0, 1);
    fuel.setNuSigmaFByGroup(2.4/16.0, 2);
    fuel.setSigmaSByGroup(71.0/360.0, 1, 1);
    fuel.setSigmaSByGroup(.02, 1, 2);
    fuel.setSigmaSByGroup(0.0, 2, 1);
    fuel.setSigmaSByGroup(11.0/15.0, 2, 2);
    fuel.setChiByGroup(1.0, 1);
    fuel.setChiByGroup(0.0, 2);
    Material* fuel_point = &fuel;

    // create moderator
    Material moderator(0, "moderator");
    moderator.setNumEnergyGroups(num_groups);
    moderator.setSigmaTByGroup(2.0/9.0, 1);
    moderator.setSigmaTByGroup(5.0/3.0, 2);
    moderator.setSigmaFByGroup(0.0, 1);
    moderator.setSigmaFByGroup(0.0, 2);
    moderator.setNuSigmaFByGroup(0.0, 1);
    moderator.setNuSigmaFByGroup(0.0, 2);
    moderator.setSigmaSByGroup(71.0/360.0, 1, 1);
    moderator.setSigmaSByGroup(.025, 1, 2);
    moderator.setSigmaSByGroup(0.0, 2, 1);
    moderator.setSigmaSByGroup(47.0/30.0, 2, 2);
    moderator.setChiByGroup(0.0, 1);
    moderator.setChiByGroup(0.0, 2);
    Material* moderator_point = &moderator;

    // create cells
    Cell root_cell(0, "root");
    root_cell.addSurface(1, x_min_point);
    root_cell.addSurface(-1, x_max_point);
    root_cell.addSurface(1, y_min_point);
    root_cell.addSurface(-1, y_max_point);
    root_cell.addSurface(1, z_min_point);
    root_cell.addSurface(-1, z_max_point);
    Cell* root_cell_point = &root_cell;
    
    Cell moderator_cell(1, "moderator");
    moderator_cell.setFill(moderator_point);
    Cell* moderator_cell_point = &moderator_cell;
    
    Cell fuel_cell(1, "fuel");
    fuel_cell.setFill(fuel_point);
    Cell* fuel_cell_point = &fuel_cell;

    // create universes
    Universe root_universe(0, "root universe");
    root_universe.addCell(root_cell_point);
    Universe moderator_universe(1, "moderator universe");
    moderator_universe.addCell(moderator_cell_point);
    Universe fuel_universe(2, "fuel universe");
    fuel_universe.addCell(fuel_cell_point);

    // create lattice
    Lattice lattice;
    lattice.setNumX(9);
    lattice.setNumY(9);
    lattice.setWidth(4.0/9.0, 4.0/9.0);
    
    // create universe array for input into lattice
    int row = 9;
    int col = 9;
    Universe** universe_array;
    universe_array = new Universe* [row];
    for (int i=0; i<row; ++i)
        universe_array[i] = new Universe [col];
    for (int i=0; i<row; ++i) {
        for (int j=0; j<col; ++j) {
            universe_array[i][j] = moderator_universe;
        }
    }
    for (int i=3; i<6; ++i) {
        for (int j=3; j<6; ++j) {
            universe_array[i][j] = fuel_universe;
        }
    }

//    std::cout << "number cells in universe " << universe_array[8][8].getNumCells()
//        << std::endl;

    // set universes in lattice
    lattice.setUniverses(1, 9, 9, universe_array);




/*
    // create geometry with surfaces
    Boundaries test_boundary;
    test_boundary.setSurface(X, MAX, &right);
    test_boundary.setSurface(X, MIN, &left);
    test_boundary.setSurface(Y, MAX, &out);
    test_boundary.setSurface(Y, MIN, &in);
    test_boundary.setSurface(Z, MAX, &top);
    test_boundary.setSurface(Z, MIN, &bottom);

    // create array with z_coordinate max and mins for 2D OpenMOC
    std::vector <double> z_bounds (2);
    z_bounds[0] = test_boundary.getSurfaceCoord(2, 0);
    z_bounds[1] = test_boundary.getSurfaceCoord(2, 1);

    

    // create mesh
    Material* point_moderator = &moderator;
    Mesh test_mesh(test_boundary, 4.0/9.0, 4.0/9.0, 4.0, point_moderator,
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
    Material* point_fuel = &fuel;
    test_mesh.fillMaterials(point_fuel, fuel_limits);

    // simulate neutron histories
    int num_neutrons = 990;
    int num_batches = 3;
    
    generateNeutronHistories(num_neutrons, test_boundary,
            test_mesh, lattice, num_batches, num_groups, z_bounds);

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
