/*
 @file      Monte_carlo.cpp
 @brief     utility functions for monte carlo neutron simulation
 @author    Luke Eure
 @date      January 12 2016
*/

#include "Monte_carlo.h"

/*
 @brief     generates and transports neutron histories, calculates the mean
            crow distance
 @param     n_histories number of neutron histories to run
 @param     mat a Material object containing information
            about the material
 @param     bounds a Boundaries object containing the limits of the
            bounding box
 @param     mesh a Mesh object containing information about the mesh
 @param     num_batches the number of batches to be tested
 @param     num_groups the number of neutron energy groups
*/
void generateNeutronHistories(int n_histories, Boundaries bounds,
        Mesh &mesh, int num_batches, int num_groups) {

    // create arrays for tallies and fissions
    std::vector <Tally> tallies(5);
    std::vector <std::vector <double> >* old_fission_bank = 
        new std::vector <std::vector <double> >;
    std::vector <std::vector <double> >* new_fission_bank = 
        new std::vector <std::vector <double> >;
    std::vector <std::vector <double> >* temp_fission_bank;
    
    bool first_round = true;

    for (int batch=1; batch <= num_batches; ++batch) {

        // clear flux data
        mesh.fluxClear();

        // assign new fission locations to old fission locations
        temp_fission_bank = old_fission_bank;
        old_fission_bank = new_fission_bank;
        new_fission_bank = temp_fission_bank;
        new_fission_bank->clear();

        // clear tallies for leaks absorptions and fissions
        tallies[LEAKS].clear();
        tallies[ABSORPTIONS].clear();
        tallies[FISSIONS].clear();

        // simulate neutron behavior
        for (int i=0; i<n_histories; ++i) {
            transportNeutron(bounds, tallies, first_round,
                    mesh, old_fission_bank, new_fission_bank, num_groups);
        }

        // give results
        std::cout << "For batch " << batch << ", k = "
            << tallies[FISSIONS].getCount()/(tallies[LEAKS].getCount() +
                    tallies[ABSORPTIONS].getCount()) << std::endl;
        first_round = false;
    }
    double mean_crow_distance = tallies[CROWS].getCount()
        / tallies[NUM_CROWS].getCount();
    std::cout << "Mean crow fly distance = " << mean_crow_distance << std::endl;

    delete old_fission_bank;
    delete new_fission_bank;
}

/*
 @brief     function that generates a neutron and measures how 
            far it travels before being absorbed.
 @details   a neutron is created in the bounding box using
            sample_location() for the first batch and sample_fission_site()
            for the rest of the batches. It moves a distance determined by
            sample_distance(). It is then either absorbed or
            scattered as determined by sample_interaction(). When
            it is absorbed, its distance from its starting point
            is appended to crow_distances. If the absorption creates
            a fission event, the number of neutrons emited is sampled.
            The location of the fission event is added to a list of fission
            events.
 @param     bounds a Boundaries object containing the limits
            of the bounding box
 @param     tallies a dictionary containing tallies of crow distances,
            leakages, absorptions, and fissions
 @param     mesh a Mesh object containing information about the mesh
 @param     old_fission_banks containing the old fission bank
 @param     new_fission_banks containing the new fission bank
 @param     num_groups the number of neutron energy groups

*/
void transportNeutron(Boundaries bounds, std::vector <Tally> &tallies,
        bool first_round, Mesh &mesh,
        std::vector <std::vector <double> >* old_fission_bank,
        std::vector <std::vector <double> >* new_fission_bank, int num_groups) {
    const double TINY_MOVE = 1e-10;
    
    // get neutron starting poinit
    std::vector <double> neutron_starting_point;
    if (first_round) {
        neutron_starting_point = sampleLocation(bounds);
    }
    else {
        neutron_starting_point = sampleFissionSite(*old_fission_bank);
    }

    // initialize neutron
    double theta = samplePolarAngle();
    double phi = sampleAzimuthalAngle();
    Neutron neutron(neutron_starting_point, theta, phi);
     
    // get mesh cell
    std::vector <double> neutron_direction;
    neutron_direction = neutron.getDirectionVector();
    std::vector <int> cell;
    cell = mesh.getCell(neutron_starting_point, neutron_direction);
    neutron.setCell(cell);

    // set neutron group
    // making mat.getChi() accept a group parameter made this much more
    // complicated
    Material* cell_mat;
    int group;
    cell_mat = mesh.getMaterial(cell);
    std::vector <double> chi(num_groups);
    for (int g=0; g<num_groups; ++g) {
        chi[g] = cell_mat->getChi(g);
    }
    group = sampleNeutronEnergyGroup(chi);
    neutron.setGroup(group);
    
    // follow neutron while it's alive
    while (neutron.alive()) {

        cell_mat = mesh.getMaterial(cell);
        group = neutron.getGroup();
        double neutron_distance;
        neutron_distance = sampleDistance(cell_mat, group);
        std::vector <double> neutron_position;
    
        // track neutron until collision or leakage
        while (neutron_distance > 0) {
            neutron_position = neutron.getPositionVector();

            // get cell boundaries
            std::vector <double> cell_mins;
            std::vector <double> cell_maxes;
            cell_mins = mesh.getCellMin(cell);
            cell_maxes = mesh.getCellMax(cell);

            // calculate distances to cell boundaries
            std::vector <std::vector <double> > distance_to_cell_edge(3);
            for (int axis=0; axis<3; ++axis) {
                distance_to_cell_edge[axis].resize(2);
                distance_to_cell_edge[axis][0] =
                    cell_mins[axis] - neutron.getPosition(axis);
                distance_to_cell_edge[axis][1] =
                    cell_maxes[axis] - neutron.getPosition(axis);
            }

            // tempd contains the current smallest r
            double tempd;
            tempd = neutron_distance;

            // create lim_bounds
            std::vector <int> cell_lim_bound;
            std::vector <int> box_lim_bound;
           
            // clear lim_bounds
            cell_lim_bound.clear();
            box_lim_bound.clear();

            // test each boundary
            double r;
            for (int axis=0; axis<3; ++axis) {
                for (int side=0; side<2; ++side) {

                    // r is variable that contains the distance
                    // along the direction vector to the boundary being tested.
                    r = distance_to_cell_edge[axis][side]
                        / neutron.getDirection(axis);
                    if (r > 0 & r < tempd) {
                        tempd = r;
                        cell_lim_bound.clear();
                        cell_lim_bound.push_back(axis*2+side);
                    }
                    else if (r == tempd) {
                        cell_lim_bound.push_back(axis*2+side);
                    }
                }
            }

            // move neutron
            neutron.move(tempd);

            // add distance to cell flux
            mesh.fluxAdd(cell, tempd, group);

            // determine boundary status
            for (int sur_side=0; sur_side <6; ++sur_side) {
                int axis = sur_side/2;
                int side = sur_side%2;

                // if sur_side is in cell_lim_bound
                if (std::find(cell_lim_bound.begin(),
                            cell_lim_bound.end(),sur_side)
                        != cell_lim_bound.end()) {
                    if (cell_mins[axis] == bounds.getSurfaceCoord(axis, side)
                            | cell_maxes[axis] == 
                            bounds.getSurfaceCoord(axis, side)) {
                        box_lim_bound.push_back(sur_side);
                    }
                }
            }

            // check boundary conditions on all hit surfaces
            for (int sur_side=0; sur_side <6; ++sur_side) {
                int axis = sur_side/2;
                int side = sur_side%2;

                // if sur_side is in box_lim_bound
                if (std::find(box_lim_bound.begin(),
                            box_lim_bound.end(),sur_side)
                        != box_lim_bound.end()) {

                    // if the neutron is reflected
                    if (bounds.getSurfaceType(axis, side) == 1) {
                        neutron.reflect(axis);

                        // place neutron on boundary to eliminate 
                        //    roundoff error
                        double bound_val;
                        bound_val = bounds.getSurfaceCoord(axis, side);
                        neutron.setPosition(axis, bound_val);
                    }

                    // if the neutron escapes
                    if (bounds.getSurfaceType(axis, side) == 0) {
                        neutron.kill();
                        neutron_distance = tempd;
                        tallies[LEAKS].add(1);
                    }
                }
            }
 
            // shorten neutron distance to collision
            neutron_distance -= tempd;

            // get new neutron cell
            if (neutron_distance > 0.0) {
                neutron_direction = neutron.getDirectionVector();
                cell = mesh.getCell(neutron_position, neutron_direction);

                // nudge neutron and find its cell
                neutron.move(TINY_MOVE);
                neutron_position = neutron.getPositionVector();
                if (mesh.positionInBounds(neutron_position)) {
                    cell = mesh.getCell(neutron_position, neutron_direction);
                }
                neutron.move(-TINY_MOVE);
                neutron.setCell(cell);
            }
        }

        // check interaction
        if (neutron.alive()) {
            cell_mat = mesh.getMaterial(cell);

            // sample what the interaction will be
            int neutron_interaction;
            neutron_interaction = sampleInteraction(cell_mat, group);

            // scattering event
            if (neutron_interaction == 0) {

                // sample scattered direction
                theta = samplePolarAngle();
                phi = sampleAzimuthalAngle();

                // sample new energy group
                int new_group;
                std::vector <double> temp_sigma_s_group;
                temp_sigma_s_group = cell_mat->getSigmaS(group);
                new_group = sampleScatteredGroup(temp_sigma_s_group, group);

                // set new group
                neutron.setGroup(new_group);

                // set new direction
                neutron.setDirection(theta, phi);
            }

            // absorption event
            else {

                // tally absorption
                tallies[ABSORPTIONS].add(1);

                // sample for fission event
                group = neutron.getGroup();
                cell = neutron.getCell();
                cell_mat = mesh.getMaterial(cell);
                neutron_position = neutron.getPositionVector();

                // fission event
                if (sampleFission(cell_mat, group) == 1) {

                    // sample number of neutrons
                    for (int i=0; i<sampleNumFission(cell_mat); ++i) {
                        new_fission_bank->push_back(neutron_position);
                        tallies[FISSIONS].add(1);
                    }
                }

                // end neutron history
                neutron.kill();
            }
        }
    }

    // tally crow distance
    double crow_distance;
    crow_distance = neutron.getDistance(neutron_starting_point);
    tallies[CROWS].add(crow_distance);
    tallies[NUM_CROWS].add(1);
}
