#!/usr/local/bin/python3
# -*- coding: utf-8 -*-


import os 
import argparse

import numpy as np


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.image import imread 

from ase import Atoms

from ase.io import read, write 
from ase.visualize import view

"""
paper figure 1. 
structure file, surface_model.con
"""

def write_surface_png(atoms, rot, povname):
    # running index for the bonds 
    bondatoms = []

    # set texttures 
    textures = ['ase3' for a in atoms]

    # set rotation 
    
    # Common kwargs for eps, png, pov
    kwargs = {
        'rotation'      : rot, # text string with rotation (default='' )
        'radii'         : .85, # float, or a list with one float per atom
        'colors'        : None,# List: one (r, g, b) tuple per atom
        'show_unit_cell': 0,   # 0, 1, or 2 to not show, show, and show all of cell
        'celllinewidth' : 0.1,  # Radius of the cylinders representing the cell
        'bondatoms'     : bondatoms, # list of tuples 
        'bondlinewidth' : 0.2, # linewidth of bond 
        'textures'      : textures, # Length of atoms list of texture names
        }
    
    # Extra kwargs only available for povray (All units in angstrom)
    kwargs.update({
        'run_povray'   : True, # Run povray or just write .pov + .ini files
        'display'      : False,# Display while rendering
        'pause'        : True, # Pause when done rendering (only if display)
        'transparent'  : True,# Transparent background
        'canvas_width' : None, # Width of canvas in pixels
        'canvas_height': 100, # Height of canvas in pixels 
        'camera_dist'  : 100.,  # Distance from camera to front atom
        'image_plane'  : None, # Distance from front atom to image plane
        'camera_type'  : 'perspective', # perspective, ultra_wide_angle
        'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
        'area_light'   : [(2., 3., 40.), # location
                          'White',       # color
                          .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
        'background'   : 'Clear',        # color
        })
       
    # Write the .pov (and .ini) file. If run_povray=False, you must run command
    # `povray filename.ini` to convert .pov file to .png
    write(povname+'.pov', atoms, **kwargs)

if __name__ == '__main__':
    # read atoms
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--atom')
    args = parser.parse_args()

    atom = Atoms(args.atom)
    write_surface_png(atom, '0x,0y,0z', args.atom)
