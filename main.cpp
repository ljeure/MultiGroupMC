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
//#include "../../OpenMOC/src/Surface.cpp"

int main() {

    // create surfaces
    MCSurface x_max(MCVACUUM, 2.0);
    MCSurface x_min(MCVACUUM, -2.0);
    MCSurface y_max(MCVACUUM, 2.0);
    MCSurface y_min(MCVACUUM, -2.0);
    MCSurface z_max(MCREFLECTIVE, 2.0);
    MCSurface z_min(MCREFLECTIVE, -2.0);

    // create openmoc surfaces and set their boundary types
    XPlane left(-2.0, 0, "left");
    XPlane right(2.1, 1, "right");
    YPlane in(-2.0, 2, "in");
    YPlane out(2.0, 3, "out");
    ZPlane bottom(-2.0, 4, "bottom");
    ZPlane top(2.0, 5, "top");
    left.setBoundaryType(VACUUM);
    right.setBoundaryType(VACUUM);
    in.setBoundaryType(VACUUM);
    out.setBoundaryType(VACUUM);
    top.setBoundaryType(REFLECTIVE);
    bottom.setBoundaryType(REFLECTIVE);

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

    // initialize lattice
    Lattice lattice;
    lattice.setNumX(9);
    lattice.setNumY(9);
    lattice.setWidth(4.0/9.0, 4.0/9.0);

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

    std::cout << std::endl;
    return 0;

    }
