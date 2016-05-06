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
 @param     flux a Flux object containing information about the flux
 @param     num_batches the number of batches to be tested
 @param     num_groups the number of neutron energy groups
*/
void generateNeutronHistories(int n_histories, Boundaries bounds,
        Flux &flux, Lattice* lattice, int num_batches, int num_groups,
        Geometry* geometry, Universe* root_universe) {

    // initialize fsid
    geometry->initializeFSRs();
    
    // create arrays for tallies and fissions
    std::vector <Tally> tallies(5);
    Fission fission_banks;
    
    bool first_round = true;

    for (int batch=1; batch <= num_batches; ++batch) {

        // clear flux data
        flux.clear();

        // assign new fission locations to old fission locations
        fission_banks.newBatch();

        // clear tallies for leaks absorptions and fissions
        tallies[LEAKS].clear();
        tallies[ABSORPTIONS].clear();
        tallies[FISSIONS].clear();

        // simulate neutron behavior
        for (int i=0; i<n_histories; ++i) {
            transportNeutron(bounds, tallies, first_round, flux, lattice,
                    &fission_banks, num_groups, i, root_universe,
                    geometry);
        }

        // give results
        double k = tallies[FISSIONS].getCount() /
            (tallies[LEAKS].getCount() + tallies[ABSORPTIONS].getCount());
        double sumStandardDev = tallies[LEAKS].getStandardDeviation(n_histories)
            + tallies[FISSIONS].getStandardDeviation(n_histories)
            + tallies[ABSORPTIONS].getStandardDeviation(n_histories);
        std::cout << "\nFor batch " << batch << ", k = " << k
            << " with standard deviation of " << sumStandardDev << std::endl;
        double fissions;
        fissions = tallies[FISSIONS].getCount();
        std::cout << "fissions: " << fissions/fissions << std::endl;
        std::cout << "leaks: " << tallies[LEAKS].getCount()/fissions
            << std::endl;
        std::cout << "absorptions: " << tallies[ABSORPTIONS].getCount()/fissions
            << std::endl;
            first_round = false;
    }
    double mean_crow_distance = tallies[CROWS].getCount()
        / tallies[NUM_CROWS].getCount();
    std::cout << "Mean crow fly distance = " << mean_crow_distance << std::endl;

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
 @param     flux a Flux object containing information about the flux
 @param     old_fission_banks containing the old fission bank
 @param     new_fission_banks containing the new fission bank
 @param     num_groups the number of neutron energy groups
 @param     neutron_num an int used for indexing neutrons
 @param     z_bounds a vector containing the minimum and maximum geometry
            boundaries along the z axis since OpenMOC only handles 2D geometries
*/
void transportNeutron(Boundaries bounds, std::vector <Tally> &tallies,
        bool first_round, Flux &flux, Lattice* lattice, Fission* fission_banks,
        int num_groups, int neutron_num, Universe* root_universe,
        Geometry* geometry) {

    const double BOUNDARY_ERROR = 1e-10;
    
    Point* neutron_position = new Point();
    
    // new way to sample neutron and set its direction
    Neutron neutron(neutron_num);
    neutron.sampleDirection();
    
    // get and set neutron starting poinit
    if (first_round)
        bounds.sampleLocation(&neutron);
    else {
        fission_banks->sampleSite(&neutron);
    }

    // get lattice cell
    Point* neutron_starting_point;
    neutron.getPositionVector(neutron_starting_point);

    std::vector <int> cell (3);
    cell[0] = lattice->getLatX(neutron_starting_point);
    cell[1] = lattice->getLatY(neutron_starting_point);
    //cell[2] = lattice->getLatZ(neutron_starting_point);
    //std::cout << "got lattice cell\n";
    
    neutron.setCell(cell);

    // set neutron group
    Material* cell_mat;
    Cell* cell_obj;
    int group;

    // get cell material
    LocalCoords* neutron_coord_position = new LocalCoords(
            neutron_starting_point->getX(), neutron_starting_point->getY(),
            neutron_starting_point->getZ());
    neutron_coord_position->setUniverse(root_universe);
    cell_obj = geometry->findCellContainingCoords(neutron_coord_position);
    cell_mat = cell_obj->getFillMaterial();
  
    std::vector <double> chi(num_groups);
    for (int g=0; g<num_groups; ++g) {
        chi[g] = cell_mat->getChiByGroup(g+1);
    }
    group = neutron.sampleNeutronEnergyGroup(chi);
    neutron.setGroup(group);
    
    // follow neutron while it's alive
    while (neutron.alive()) {

// stuff from segmetize --------------------------------------------------------

        // use a LocalCoords for the start and end of each segment
        double neutron_distance;
            
        //sample a distance to travel in the material
        neutron_distance = 
            -log(neutron.arand()) / cell_mat->getSigmaTByGroup(group+1);
    
        // track neutron until collision or escape
        while (neutron_distance > 0.000001) {

            double x0 = neutron_position->getX();
            double y0 = neutron_position->getY();
            double z0 = neutron_position->getZ();
            double phi =
                atan(neutron.getDirection(1)/neutron.getDirection(0));
            FP_PRECISION length;
            Material* material;
            int fsr_id;

            LocalCoords start(x0, y0, z0);
            LocalCoords end(x0, y0, z0);
            start.setUniverse(root_universe);
            end.setUniverse(root_universe);
            start.setPhi(phi);
            end.setPhi(phi);
    
            // find the Cell containing the Track starting Point 
            Cell* curr = geometry->findFirstCell(&end);
            Cell* prev;

            // if starting Point was outside the bounds of the Geometry /
            if (curr == NULL)
                log_printf(ERROR,
                        "Could not find a material-filled Cell containing"
                        " the start Point of this Track: %s");

            // While the end of the segment's LocalCoords is still
            // within the Geometry move it to the next Cell, create
            // a new segment, and add it to the Geometry
            while (curr != NULL) {

                end.copyCoords(&start);
                end.setPhi(phi);

                // Find the next Cell along the Track's trajectory
                prev = curr;
                curr = geometry->findNextCell(&end);

                // Checks that segment does not have the same
                // start and end Points 
                if (start.getX() == end.getX()
                        && start.getY() == end.getY())
                    log_printf(ERROR,
                            "Created segment with same start and end point:"
                            " x = %f, y = %f",start.getX(), start.getY());

                // Find the segment length, Material and
                // FSR ID of the last cell
                length =
                    FP_PRECISION(end.getPoint()->distanceToPoint(
                                start.getPoint()));
                //std::cout << "length " << length << std::endl;
                fsr_id = geometry->findFSRId(&start);
    
                            std::cout << "end " << end.getX() << " "
                                << end.getY() << " "
                                << end.getZ() << std::endl;
                            std::cout << "start " << start.getX() << " "
                                << start.getY() << " "
                                << start.getZ() << std::endl;
                // if neutron's path doesn't end in this cell
                if (length < neutron_distance) {

                    // add distance travelled to flux, shorten distance and
                    // move neutron
                    flux.add(neutron.getGroup(), fsr_id, length);
                    neutron_distance -= length;
                    neutron.setPosition(0, end.getX());
                    neutron.setPosition(1, end.getY());
                    neutron.setPosition(2, end.getZ());
                }

                // if the neutron's path ends in this cell
                else {
                    flux.add(neutron.getGroup(), fsr_id, length);
                    neutron.move(neutron_distance);
                    for (int i=0; i<2; i++) {
                        if (neutron.getPosition(i) > 2.0001
                                or neutron.getPosition(i) < -2.0001) {
/*
                            std::cout << "end " << end.getX() << " "
                                << end.getY() << " "
                                << end.getZ() << std::endl;
                            std::cout << "start " << start.getX() << " "
                                << start.getY() << " "
                                << start.getZ() << std::endl;
                            std::cout << "neutron distance "
                                << neutron_distance << std::endl;
                            std::cout << "length "
                                << length << std::endl;
                            neutron.move(-neutron_distance);
                            std::cout
                                << "----------------------premove position "
                                << neutron.getPosition(0) <<
                                " " << neutron.getPosition(1) <<
                                " " << neutron.getPosition(2) << std::endl;
                            neutron.move(neutron_distance);
                            std::cout
                                << "------------------------------position "
                                << neutron.getPosition(0) <<
                                " " << neutron.getPosition(1) <<
                                " " << neutron.getPosition(2) << std::endl;
   */
                         }
                    }
                    neutron_distance = 0;
                    break;
                }
            } // while curr != null

            // find out which boundary the neutron is on
            std::vector <int> box_lim_bound;
            box_lim_bound.clear();
            if (std::abs(neutron.getPosition(0)
                        - bounds.getSurfaceCoord(0, 0)) < BOUNDARY_ERROR) {
                // std::cout << "hit boundary 0\n";
                box_lim_bound.push_back(0);
            }
            if (std::abs(neutron.getPosition(0)
                        - bounds.getSurfaceCoord(0, 1)) < BOUNDARY_ERROR) { 
               // std::cout << "hit boundary 1\n";
                box_lim_bound.push_back(1);
            }
            if (std::abs(neutron.getPosition(1)
                        - bounds.getSurfaceCoord(1, 0)) < BOUNDARY_ERROR) {
               // std::cout << "hit boundary 2\n";
                box_lim_bound.push_back(2);
            }
            if (std::abs(neutron.getPosition(1)
                        - bounds.getSurfaceCoord(1, 1)) < BOUNDARY_ERROR) {
               // std::cout << "hit boundary 3\n";
                box_lim_bound.push_back(3);
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

                        // place neutron on boundary to 
                        // eliminate roundoff error
                        double bound_val;
                        bound_val = bounds.getSurfaceCoord(axis, side);
                        neutron.setPosition(axis, bound_val);
                    }

                    // if the neutron escapes
                    if (bounds.getSurfaceType(axis, side) == 0) {
                        neutron.kill();
                        neutron_distance = 0;
                        tallies[LEAKS] += 1;
                        break;
                    }
                }
            } // check boundary conditions
        } // while distance > 0

        // this is where to go when distance is up

//-----------------------------------------------------------------------------
/*

        // track neutron until collision or leakage
        while (neutron_distance > BOUNDARY_ERROR) {
            neutron.getPositionVector(neutron_position);

            // get cell boundaries
            std::vector <double> cell_mins (3);
            std::vector <double> cell_maxes (3);
            cell_mins[0] = lattice->getMinX() + cell[0] * lattice->getWidthX();
            cell_mins[1] = lattice->getMinY() + cell[1] * lattice->getWidthY();
            cell_mins[2] = lattice->getMinZ() + cell[2] * lattice->getWidthZ();
            std::cout << "cell mins in z " << cell_mins[2] << std::endl;
            cell_maxes[0] = lattice->getMinX() + 
                (cell[0] + 1) * lattice->getWidthX();
            cell_maxes[1] = lattice->getMinY() + 
                (cell[1] + 1) * lattice->getWidthY();
            cell_maxes[2] = lattice->getMinZ() + 
                (cell[2] + 1) * lattice->getWidthZ();

            // calculate distances to cell boundaries
            std::vector <std::vector <double> > distance_to_cell_edge(3);
            for (int axis=0; axis<3; ++axis) {
                distance_to_cell_edge[axis].resize(2);
                distance_to_cell_edge[axis][0] =
                    cell_mins[axis] - neutron.getPosition(axis);
                distance_to_cell_edge[axis][1] =
                    cell_maxes[axis] - neutron.getPosition(axis);
            }

            // create lim_bounds
            std::vector <int> cell_lim_bound;
            std::vector <int> box_lim_bound;
           
            // clear lim_bounds
            cell_lim_bound.clear();
            box_lim_bound.clear();

            // tempd contains the current smallest r
            double tempd;
            tempd = neutron_distance;

            // test each boundary
            double r;
            for (int axis=0; axis<3; ++axis) {
                for (int side=0; side<2; ++side) {

                    // r is variable that contains the distance
                    // along the direction vector to the boundary being tested.
                    r = distance_to_cell_edge[axis][side]
                        / neutron.getDirection(axis);

                    if (r > BOUNDARY_ERROR & r < tempd) {
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

            // add distance to fsr flux
            neutron_coord_position->setX((cell_mins[0] + cell_maxes[0])/2.);
            neutron_coord_position->setY((cell_mins[1] + cell_maxes[1])/2.);
            neutron_coord_position->setZ((cell_mins[2] + cell_maxes[2])/2.);
            neutron_coord_position->setUniverse(root_universe);
            geometry->findFirstCell(neutron_coord_position);
            int fsr = geometry->findFSRId(neutron_coord_position);
            flux.add(neutron.getGroup(), fsr, tempd);

            // shorten neutron distance to collision
            neutron_distance -= tempd;

            // determine potential geometry boundaries
            for (int sur_side=0; sur_side <4; ++sur_side) {
                int axis = sur_side/2;
                int side = sur_side%2;

                // if sur_side is in cell_lim_bound
                if (std::find(cell_lim_bound.begin(),
                            cell_lim_bound.end(),sur_side)
                        != cell_lim_bound.end()) {
                    
                    if (neutron.getCell()[axis] == 0 && side ==0)
                        box_lim_bound.push_back(sur_side);

                    // if the neutron is on the max side of the cell and the
                    // cell is the highest cell in the geometry
                    if (axis == 0) {
                        if (neutron.getCell()[axis] == lattice->getNumX() - 1
                                && side == 1)
                            box_lim_bound.push_back(sur_side);
                    }
                    if (axis == 1) {
                        if (neutron.getCell()[axis] == lattice->getNumY() - 1
                                && side == 1)
                            box_lim_bound.push_back(sur_side);
                    }
                    if (axis == 2) {
                        if (neutron.getCell()[axis] == lattice->getNumZ() - 1
                                && side == 1)
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

                        // place neutron on boundary to eliminate roundoff error
                        double bound_val;
                        bound_val = bounds.getSurfaceCoord(axis, side);
                        neutron.setPosition(axis, bound_val);
                    }

                    // if the neutron escapes
                    if (bounds.getSurfaceType(axis, side) == 0) {
                        neutron.kill();
                        neutron_distance = 0;
                        tallies[LEAKS] += 1;
                        break;
                    }
                }
            }
 
            // get new neutron cell
            if (neutron_distance > 0) {
                cell[0] = lattice->getLatX(neutron_position);
                cell[1] = lattice->getLatY(neutron_position);
                //cell[2] = lattice->getLatZ(neutron_position);
                
                // nudge neutron and find its cell
                neutron.move(1e-3);
                neutron.getPositionVector(neutron_position);

                // withinBounds seems to return true if the point is without
                if (lattice->withinBounds(neutron_position)) {
                    cell[0] = lattice->getLatX(neutron_position);
                    cell[1] = lattice->getLatY(neutron_position);
                    //cell[2] = lattice->getLatZ(neutron_position);
                }
                neutron.move(-1e-3);
  
                neutron.setCell(cell);
            }
        }
    */

        // check interaction
        if (neutron.alive()) {
            neutron_coord_position->setX(neutron_position->getX());
            neutron_coord_position->setY(neutron_position->getY());
            neutron_coord_position->setZ(neutron_position->getZ());
            neutron_coord_position->setUniverse(root_universe);
            cell_obj =
                geometry->findCellContainingCoords(neutron_coord_position);
            cell_mat = cell_obj->getFillMaterial();

            // calculate sigma_s for a group in order to sample an interaction
            std::vector <double> sigma_s_group;
            double sum_sigma_s_group=0;
            for (int g=1; g<=cell_mat->getNumEnergyGroups(); ++g) {
                sigma_s_group.push_back(cell_mat->getSigmaSByGroup(group+1, g));
                sum_sigma_s_group += cell_mat->getSigmaSByGroup(group+1, g);
            }

            // calculate sigma_a in order to sample an interaction
            double sigma_a;
            sigma_a =
                cell_mat->getSigmaTByGroup(group+1) - sum_sigma_s_group;

            // sample an interaction
            int neutron_interaction =
                (int) (neutron.arand() < (sigma_a
                            / cell_mat->getSigmaTByGroup(group+1)));
            
            // scattering event
            if (neutron_interaction == 0) {

                // sample scattered direction
                neutron.sampleDirection();

                // sample new energy group
                int new_group = neutron.sampleScatteredGroup(sigma_s_group,
                    group);

                // set new group
                neutron.setGroup(new_group);
            }

            // absorption event
            else {

                // tally absorption
                tallies[ABSORPTIONS] += 1;

                // sample for fission event
                group = neutron.getGroup();
                cell = neutron.getCell();
                neutron_coord_position->setX(neutron_position->getX());
                neutron_coord_position->setY(neutron_position->getY());
                neutron_coord_position->setZ(neutron_position->getZ());
                neutron_coord_position->setUniverse(root_universe);
                cell_obj =
                    geometry->findCellContainingCoords(neutron_coord_position);
                cell_mat = cell_obj->getFillMaterial();
                neutron.getPositionVector(neutron_position);

                // sample whether fission event occurs
                int fission_occurs =
                    neutron.arand() < cell_mat->getSigmaFByGroup(group+1)
                    / sigma_a;

                // fission event
                if (fission_occurs == 1) {

                    // sample number of neutrons released during fission
                    double nu = cell_mat->getNuSigmaFByGroup(1)
                        / cell_mat->getSigmaFByGroup(1);
                    int lower = (int) nu;
                    int add = (int) (neutron.arand() < nu -lower);
                    int num_neutrons = lower + add;
                    for (int i=0; i<num_neutrons; ++i) {
                        fission_banks->add(neutron_position);
                        tallies[FISSIONS] += 1;
                    }
                }

                // end neutron history
                neutron.kill();
            } // absorption event
        } // if neutron alive
    } // while neutron alive
    
    delete neutron_coord_position;

    // tally crow distance
    double crow_distance;
    crow_distance = neutron.getDistance(neutron_starting_point);
    tallies[CROWS] += crow_distance;
    tallies[NUM_CROWS] += 1;
}





/*

cal coords set phi should be the azimuthal angle
429 
430     basically use what's in geometry.segmetize
431     copy and paste from Use a local coorsd fro the statrt and end
432     to create a new track segment
433 
434     create neutron
435     curr = find first cell
436     while distance < total distance left to travel and neutron alive
437         while curr != null
438         LocaCoords end;
439         find next cell
440 
441         check boundaries
442             reflect/kill
443 
444     look at cpusolver.h and .cpp when making plot.
445     Monte carlo solver will be subclass of solver
446 

*/
