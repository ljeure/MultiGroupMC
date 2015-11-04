'''
 @file      tester.py
 @brief     A teting program for debugging
 @author    Luke Eure
 @date      November 1, 2015
'''

#import libraries
import boundaries
import coords
import distributions
import material
import monte_carlo
import plotter
import tally
import fission
import mesh

bounds =  boundaries.Boundaries(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

test_mesh = mesh.Mesh(bounds, .1)

test_material = material.Material(2.0, 1.0, 2.4, .5)

monte_carlo.generate_neutron_histories(1000, test_material, bounds, test_mesh, 1)
