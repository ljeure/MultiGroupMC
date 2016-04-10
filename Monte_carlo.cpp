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
        Mesh &mesh, Lattice &lattice, int num_batches, int num_groups,
        Geometry* geometry, std::vector <double> &z_bounds) {

    // initialize fsid
    geometry->initializeFSRs();
    
    // test fsrid
    LocalCoords* test = new LocalCoords(.1,.1,.1);
    Universe* fuel_uni = geometry->getAllUniverses()[1];
    test->setUniverse(fuel_uni);
    std::cout << "got here\n";
    Cell* fuel_cell = fuel_uni->getCell(1);
    test->setCell(fuel_cell);
    std::cout << "fsr test " << geometry->findFSRId(test) << std::endl;

    // add FSR points to geometry
    Point sampleFSR;
    for (int i=0; i<lattice.getNumX(); ++i) {
        double xPos = lattice.getMinX() + lattice.getWidthX()*(i + 1/2);
        sampleFSR.setX(xPos);
        
        for (int j=0; j<lattice.getNumX(); ++j) {
            double yPos = lattice.getMinY() + lattice.getWidthY()*(j + 1/2);
            sampleFSR.setY(yPos);
        
            for (int k=0; k<lattice.getNumZ(); ++k) {
                double zPos = lattice.getMinZ() + lattice.getWidthZ()*(k + 1/2);
                sampleFSR.setZ(zPos);
                std::cout << "k " << k << std::endl;
                LocalCoords* coordFSR = new LocalCoords(
                        sampleFSR.getX(), sampleFSR.getY(), sampleFSR.getX());
                
                std::cout << "fsrid " << i << " " << j << " " << k
                    << ": " << /*geometry->findFSRId(coordFSR) <<*/ std::endl;

            }
        }
    }



    // create arrays for tallies and fissions
    std::vector <Tally> tallies(5);
    Fission fission_banks;
    
    bool first_round = true;

    for (int batch=1; batch <= num_batches; ++batch) {

        // clear flux data
        mesh.fluxClear();

        // assign new fission locations to old fission locations
        fission_banks.newBatch();

        // clear tallies for leaks absorptions and fissions
        tallies[LEAKS].clear();
        tallies[ABSORPTIONS].clear();
        tallies[FISSIONS].clear();

        // simulate neutron behavior
        for (int i=0; i<n_histories; ++i) {
            transportNeutron(bounds, tallies, first_round, mesh, lattice,
                    &fission_banks, num_groups, i, z_bounds);
        }

        // give results
        double k = tallies[FISSIONS].getCount() /
            (tallies[LEAKS].getCount() + tallies[ABSORPTIONS].getCount());
        double sumStandardDev = tallies[LEAKS].getStandardDeviation(n_histories)
            + tallies[FISSIONS].getStandardDeviation(n_histories)
            + tallies[ABSORPTIONS].getStandardDeviation(n_histories);
        std::cout << "For batch " << batch << ", k = " << k
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
 @param     mesh a Mesh object containing information about the mesh
 @param     old_fission_banks containing the old fission bank
 @param     new_fission_banks containing the new fission bank
 @param     num_groups the number of neutron energy groups
 @param     neutron_num an int used for indexing neutrons
 @param     z_bounds a vector containing the minimum and maximum geometry
            boundaries along the z axis since OpenMOC only handles 2D geometries
*/
void transportNeutron(Boundaries bounds, std::vector <Tally> &tallies,
        bool first_round, Mesh &mesh, Lattice &lattice, Fission* fission_banks,
        int num_groups, int neutron_num, std::vector <double> &z_bounds) {

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

    // get mesh cell
    Point* neutron_starting_point;
    neutron.getPositionVector(neutron_starting_point);

    std::vector <int> cell (3);
//    std::cout << "starting point x " << neutron_starting_point->getX() << std::endl;
    cell[0] = lattice.getLatX(neutron_starting_point);
    cell[1] = lattice.getLatY(neutron_starting_point);
    cell[2] = lattice.getLatZ(neutron_starting_point);
    
    neutron.setCell(cell);

    // set neutron group
    Material* cell_mat;
    int group;
    cell_mat = mesh.getMaterial(cell);
    std::vector <double> chi(num_groups);
    for (int g=0; g<num_groups; ++g) {
        chi[g] = cell_mat->getChiByGroup(g+1);
    }
    group = neutron.sampleNeutronEnergyGroup(chi);
    neutron.setGroup(group);
    
    // follow neutron while it's alive
    while (neutron.alive()) {

        cell_mat = mesh.getMaterial(cell);
        group = neutron.getGroup();
        double neutron_distance;

        //sample a distance to travel in the material
        neutron_distance = 
            -log(neutron.arand()) / cell_mat->getSigmaTByGroup(group+1);
    
        // track neutron until collision or leakage
        while (neutron_distance > BOUNDARY_ERROR) {
            neutron.getPositionVector(neutron_position);

            // get cell boundaries
            std::vector <double> cell_mins (3);
            std::vector <double> cell_maxes (3);
            cell_mins[0] = lattice.getMinX() + cell[0] * lattice.getWidthX();
            cell_mins[1] = lattice.getMinY() + cell[1] * lattice.getWidthY();
            cell_mins[2] = z_bounds[0];
            cell_maxes[0] = lattice.getMinX() + 
                (cell[0] + 1) * lattice.getWidthX();
            cell_maxes[1] = lattice.getMinY() + 
                (cell[1] + 1) * lattice.getWidthY();
            cell_maxes[2] = z_bounds[1];

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

     //       std::cout << "position before moving " << neutron.getPosition(0)
      //          << std::endl;
            // move neutron
            neutron.move(tempd);
 //           std::cout << "position after moving " << neutron.getPosition(0)
   //             << std::endl;

            // add distance to cell flux
            mesh.fluxAdd(cell, tempd, group);

            // shorten neutron distance to collision
            neutron_distance -= tempd;

            // determine potential geometry boundaries
            for (int sur_side=0; sur_side <6; ++sur_side) {
                int axis = sur_side/2;
                int side = sur_side%2;

                // if sur_side is in cell_lim_bound
                if (std::find(cell_lim_bound.begin(),
                            cell_lim_bound.end(),sur_side)
                        != cell_lim_bound.end()) {
                    
                    if (neutron.getCell()[axis] == 0 && side ==0)
                        box_lim_bound.push_back(sur_side);

                    // change to lattice when it gets 3D capabilities
                    if (neutron.getCell()[axis] == mesh.getNumCells(axis) - 1
                            && side == 1)
                        box_lim_bound.push_back(sur_side);
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
  /*              std::cout << "neutron position " << neutron_position->getX()
                    << std::endl;
    */
                    cell[0] = lattice.getLatX(neutron_position);
                cell[1] = lattice.getLatY(neutron_position);
                cell[2] = lattice.getLatZ(neutron_position);
                
                // nudge neutron and find its cell
                neutron.move(1e-3);
                neutron.getPositionVector(neutron_position);
                /*std::cout << "neutron position after nudge "
                    << neutron_position->getX() << std::endl;
*/
                // withinBounds seems to return true if the point is without
                if (lattice.withinBounds(neutron_position)) {
                    cell[0] = lattice.getLatX(neutron_position);
                    cell[1] = lattice.getLatY(neutron_position);
                    cell[2] = lattice.getLatZ(neutron_position);
                }
                neutron.move(-1e-3);
  
                /*std::cout << "neutron cell currently "
                    << neutron.getCell()[0] <<std::endl;
                */
                    neutron.setCell(cell);
                //std::cout << "neutron cell set to " << cell[0] <<std::endl;
            }
        }

        // check interaction
        if (neutron.alive()) {
            cell_mat = mesh.getMaterial(cell);

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
                cell_mat = mesh.getMaterial(cell);
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
            }
        }
    }

    // tally crow distance
    double crow_distance;
    crow_distance = neutron.getDistance(neutron_starting_point);
    tallies[CROWS] += crow_distance;
    tallies[NUM_CROWS] += 1;
}
