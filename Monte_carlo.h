/* 
 @file      Monte_carlo.h
 @brief     header file for monte_carlo.cpp
 @author    Luke Eure
 @date      January 12 2016
*/

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "../../OpenMOC/src/Point.h"
#include "../../OpenMOC/src/Universe.h"
#include "../../OpenMOC/src/Material.h"
#include "../../OpenMOC/src/Geometry.h"
#include "../../OpenMOC/src/LocalCoords.h"

#include "Enumerations.h"
#include "Tally.h"
#include "Mesh.h"
#include "Neutron.h"
#include "Fission.h"

enum tally_names {CROWS, NUM_CROWS, LEAKS, ABSORPTIONS, FISSIONS};
enum fission_bank_names {OLD, NEW};

void generateNeutronHistories(int n_histories, Boundaries bounds,
        Mesh &mesh, Lattice* lattice, int num_batches, int num_groups,
        Geometry* geometry, Universe* root_universe);

void transportNeutron(Boundaries bounds, std::vector <Tally> &tallies,
        bool first_round, Mesh &mesh, Lattice* lattice, Fission* fission_banks,
        int num_groups, int neutron_num, Universe* root_universe,
        Geometry* geometry);

#endif
