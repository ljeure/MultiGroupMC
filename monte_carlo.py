'''
 @file      monte_carlo.py
 @brief     Utility functions for monte carlo neutron simulations
 @details   This file contains functions to be used in a Monte Carlo
            simulation. 
            Later, these funcitons will be re-organized and implemented as 
            objects.
 @author    Geoffrey Gunow, Luke Eure
 @date      September 30, 2015
'''

# Import libraries
from math import *
from distributions import *
from plotter import *
import copy
import coords
import tally
import fission

'''
 @brief     Function that generates a neutron and measures how 
            far it travels before being absorbed.
 @details   A neutron is created in the bounding box using
            sample_location(). It moves a distance determined by
            sample_distance(). It is then either absorbed or
            scattered as determined by sample_interaction(). When
            it is absorbed, its distance from its starting point
            is appended to crow_distances.
 @param	    mat an instance of the Material class containing information
            about the material
 @param     bounds an instance of the boundaries class containing the limits
            of the bounding box
 @param     tallies a dictionary containing tallies of crow distances,
            leakages, absorptions, and fissions
 @param     fission_banks a dictionary containing the old and new fission banks
 @param     first_round a boolean telling whenther or not this is the first
            batch to be tested
'''
def transport_neutron(mat, bounds, tallies, fission_banks, first_round):
    
    '''
    for i in xrange(10):
    
        old_fission_locations.clear()
        for num in xrange(new_fission_locations.length):
            old_fission_locations.add(new_fission_locations.next())
        new_fission_locations.clear()
    '''

    # first round
    if first_round:
        neutron_starting_point = sample_location(bounds)
        
    # every other round
    else:
        neutron_starting_point = sample_fission_site(fission_banks["old"])
        
    neutron_lost = False
    neutron_reflected = False
    neutron_position = coords.Coords(neutron_starting_point.x,
            neutron_starting_point.y, neutron_starting_point.z) 
    neutron_interaction = 0
    while neutron_interaction != 1:
        neutron_distance = sample_distance(mat)
        theta = sample_polar_angle()
        phi = sample_azimuthal_angle()
    
        # neutron movement unit vector
        neutron_movement_unit_vector = coords.Coords(sin(theta)*cos(phi),
                   sin(theta)*sin(phi), cos(theta))

        while neutron_distance > 0: 
            xmaxd = bounds.get_x_max() - neutron_position.x
            xmind = bounds.get_x_min() - neutron_position.x
            ymaxd = bounds.get_y_max() - neutron_position.y
            ymind = bounds.get_y_min() - neutron_position.y
            zmaxd = bounds.get_z_max() - neutron_position.z
            zmind = bounds.get_z_min() - neutron_position.z
               
            # lim_bound contains the string of the limiting boundary
            lim_bound = list()
                
            # tempd contains the current smallest r
            tempd = neutron_distance
                
            # r is a variable that contains the distance along the
            # direction vector to the boundary being tested 
            #test xmax
            r = xmaxd/(sin(theta)*cos(phi))
            if (0 <= r and r < neutron_distance and not r == -0.0):
                tempd = r
                lim_bound.append('x_max')
                
            #test xmin
            r = xmind/(sin(theta)*cos(phi))
            if (0 <= r and r < neutron_distance and \
                    r < tempd and not r == -0.0):
                tempd = r
                lim_bound = []
                lim_bound.append('x_min')
            elif r == tempd:
                lim_bound.append('x_min')
            
            #test ymax
            r = ymaxd/(sin(theta)*sin(phi))
            if (0 <= r and r < neutron_distance and \
                    r < tempd and not r == -0.0):
                tempd = r
                lim_bound = []
                lim_bound.append('y_max')
            elif r == tempd:
                lim_bound.append('y_max')               
                
            #test ymin
            r = ymind/(sin(theta)*sin(phi))
            if (0 <= r and r < neutron_distance and \
                    r < tempd and not r == -0.0):
                tempd = r
                lim_bound = []
                lim_bound.append('y_min')
            elif r == tempd:
                lim_bound.append('y_min')               
            
            #test zmax
            r = zmaxd/cos(theta)
            if (0 <= r and r < neutron_distance and \
                    r < tempd and not r == -0.0):
                tempd = r
                lim_bound = []
                lim_bound.append('z_max')
            elif r == tempd:
                lim_bound.append('z_max')              
                
            #test zmin
            r = zmind/cos(theta)
            if (0 <= r and r < neutron_distance and \
                    r < tempd and not r == -0.0):
                tempd = r
                lim_bound = []
                lim_bound.append('z_min')
            elif r == tempd:
                lim_bound.append('z_min')               
                
            # if the neutron stays within the bounding box
            if lim_bound == []:
                neutron_distance = 0
            neutron_position.sadd(neutron_movement_unit_vector * tempd)

            for surface in lim_bound:  
                # if the neutron is reflected
                if (bounds.get_surface_type(surface) == 1):
                    if surface == 'x_max' or surface == 'x_min':
                        neutron_movement_unit_vector.smul([-1, 1, 1])
                    if surface == 'y_max' or surface == 'y_min':
                        neutron_movement_unit_vector.smul([1, -1, 1])
                    if surface == 'z_max' or surface == 'z_min':
                        neutron_movement_unit_vector.smul([1, 1, -1])
                    neutron_reflected = True
    
                # if the neutron escapes
                elif (bounds.get_surface_type(surface) == 0):
                    neutron_lost = True
                    neutron_interaction = 1
                    neutron_distance = 0
                    tallies["leaks"].increment(1)
            if neutron_reflected:
                neutron_distance -= tempd
                
        if not neutron_lost:
            neutron_interaction = sample_interaction(mat)


    # neutron_position turns into a vector pointing from the neutron's origin
    # to its absorption point
    if not neutron_lost:
        tallies["absorptions"].increment(1)
        if sample_fission(mat) == 1:
            for j in xrange(sample_num_fission(mat)):
                fission_banks["new"].add([neutron_position.x,
                        neutron_position.y, neutron_position.z])
                tallies["fissions"].increment(1) 
        neutron_position.ssub([neutron_starting_point.x,
            neutron_starting_point.y, neutron_starting_point.z])
        tallies["crows"].increment(neutron_position.getDistance())
        tallies["num_crows"].increment(1)
                            
'''
 @brief     Generates and transports neutron histories, calculates the mean
            crow distance
 @param     n_histories number of neutron histories to run
 @param     mat an instance of the Material class containing information
            about the material
 @param     bounds an instance of the boundaries class containing the limits
            of the bounding box
'''
def generate_neutron_histories(n_histories, mat, bounds):
    
    crow_distances = tally.Tally()
    num_crow_distances = tally.Tally()
    leakage_tally = tally.Tally()
    absorption_tally = tally.Tally()
    fission_tally = tally.Tally()
    
    tallies = {"crows": crow_distances, "num_crows": num_crow_distances,
            "leaks": leakage_tally, "absorptions": absorption_tally,
            "fissions": fission_tally}

    old_fission_locations = fission.Fission()
    new_fission_locations = fission.Fission()
    fission_banks = {"old": old_fission_locations,
            "new": new_fission_locations}
    
    first_round = True
    for batch in xrange(10):
        
        # assign the new fission locations to the old fission list
        old_fission_locations.clear()
        for num in xrange(new_fission_locations.length):
            old_fission_locations.add(new_fission_locations.next())
        new_fission_locations.clear()
        
        # clear the tallies for leaks absorptions and fissions
        tallies["leaks"].clear()
        tallies["absorptions"].clear()
        tallies["fissions"].clear()

        # simulate the neutron behavior
        for i in xrange(n_histories):
            transport_neutron(mat, bounds, tallies, fission_banks, first_round)
        
        
        print "For batch ", batch + 1, ", k = ", \
                tallies["fissions"].amt/(tallies["leaks"].amt + \
                tallies["absorptions"].amt)
        first_round = False
    
    mean_crow_distance = crow_distances.amt / num_crow_distances.amt
    print "Mean crow fly distance = ", mean_crow_distance

'''
 @brief     Function that determines what a neutron's new position
            is, taking into account the surface
 @details   A neutron is created in the bounding box using
            sample_location(). It moves a distance determined by
            sample_distance(). It is then either absorbed or
            scattered as determined by sample_interaction(). When
            it is absorbed, its distance from its starting point
            is appended to crow_distances.
 @param	    mat an instance of the Material class containing information
            about the material
 @param     bounds an instance of the boundaries class containing the limits
            of the bounding box
 @param     crow_distances a list containg the distances
            from each neutron's starting point to the point that
            it is absorbed.





nu = 2.4
nu is average number of neutrons released in a fission event


'''
