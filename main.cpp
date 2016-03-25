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

int main() {

    // create surfaces
    MCSurface x_max(MCVACUUM, 2.0);
    MCSurface x_min(MCVACUUM, -2.0);
    MCSurface y_max(MCVACUUM, 2.0);
    MCSurface y_min(MCVACUUM, -2.0);
    MCSurface z_max(MCREFLECTIVE, 2.0);
    MCSurface z_min(MCREFLECTIVE, -2.0);

    // create geometry with surfaces
    Boundaries test_boundary;
    test_boundary.setSurface(X, MAX, &x_max);
    test_boundary.setSurface(X, MIN, &x_min);
    test_boundary.setSurface(Y, MAX, &y_max);
    test_boundary.setSurface(Y, MIN, &y_min);
    test_boundary.setSurface(Z, MAX, &z_max);
    test_boundary.setSurface(Z, MIN, &z_min);

    // create OpenMOC materials
    int num_groups = 2;
    Material uranium(1, "uranium");
    uranium.setNumEnergyGroups(2);
    uranium.setSigmaTByGroup(2.0/9.0, 1);
    uranium.setSigmaTByGroup(5.0/6.0, 2);
    uranium.setSigmaFByGroup(1.0/480.0, 1);
    uranium.setSigmaFByGroup(1.0/16.0, 2);
    uranium.setNuSigmaFByGroup(2.4/480.0, 1);
    uranium.setNuSigmaFByGroup(2.4/16.0, 2);
    uranium.setSigmaSByGroup(71.0/360.0, 1, 1);
    uranium.setSigmaSByGroup(.02, 1, 2);
    uranium.setSigmaSByGroup(0.0, 2, 1);
    uranium.setSigmaSByGroup(11.0/15.0, 2, 2);
    uranium.setChiByGroup(1.0, 1);
    uranium.setChiByGroup(0.0, 2);

    Material moderator(0, "moderator");
    moderator.setNumEnergyGroups(2);
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
    Material* point_uranium = &uranium;
    test_mesh.fillMaterials(point_uranium, fuel_limits);

    // initialize lattice
    Lattice lattice;
    lattice.setNumX(9);
    lattice.setNumY(9);
    lattice.setWidth(4.0/9.0, 4.0/9.0);

    // simulate neutron histories
    int num_neutrons = 1000;
    int num_batches = 3;
    
    generateNeutronHistories(num_neutrons, test_boundary,
            test_mesh, lattice, num_batches, num_groups);

    // plot neutron flux
    std::vector <std::vector <std::vector <std::vector <double> > > > flux =
       test_mesh.getFlux();
    printFluxToFile(flux);

    // run python script to get flux plots
    system("python Flux_parser.py");

    std::cout << std::endl;
    return 0;

    }
