#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

import os 
import warnings

import argparse

import json

import numpy as np
norm = np.linalg.norm

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.image import imread 

from ase import Atoms
from ase.visualize import view
from ase.io import read, write 
from ase.io.pov import get_bondpairs ,set_high_bondorder_pairs

from ase.symbols import Symbols, symbols2numbers
s2n = symbols2numbers

description="""
How to plot atoms with PovRay.
"""


def write_pov_png(atoms, povname, params):
    # running index for the bonds 
    natoms = len(atoms)

    # check 
    #atom_list = []
    #cu_pos = np.array([7.8,7.0,6.4])
    #for i, atom in enumerate(atoms):
    #    if norm(atom.position-cu_pos) < 9.0:
    #        if atom.symbol != 'C' and atom.symbol != 'H':
    #            atom_list.append(i)
    #    if norm(atom.position-cu_pos) < 4.5:
    #        if atom.symbol == 'C' or atom.symbol == 'H':
    #            atom_list.append(i)
    #atoms = atoms[atom_list]


    # check custom radii
    custom_radii = []

    radii = params.pop('atom_radii', None)
    if radii:
        for sym in atoms.get_chemical_symbols():
            if sym in radii.keys():
                custom_radii.append(radii[sym])
            else:
                warnings.warn('No %s in custom radii and use default 1.0 radius.' %sym, RuntimeWarning)
                custom_radii.append(1.0)
    else:
        custom_radii = [1.0 for at in atoms]

    print(custom_radii)

    # colors 
    custom_colors = []

    colors = params.pop('colors', None)
    if colors:
        for sym in atoms.get_chemical_symbols():
            if sym in colors.keys():
                custom_colors.append(colors[sym])
            else:
                warnings.warn('No %s in custom color and use default.' %sym, RuntimeWarning)
                custom_radii.append(None)
    else:
        custom_colors = None

    # rotation
    custom_rot = '{}x,{}y,{}z'.format(*params['rotation']) # using ag: 'view -> rotate'

    # bond pairs 
    bond_radius = params.pop('bond_times', 1.0)
    bondpairs = get_bondpairs(atoms, radius=bond_radius)
    #print(bondpairs)

    high_order_bonds = params.pop('high_order_bonds', None)
    if high_order_bonds:
        high_bondorder_pairs = {}
        for key, value in high_order_bonds.items():
            new_key = key.split('-')
            high_bondorder_pairs[(int(new_key[0]),int(new_key[1]))] = value
        #print(high_bondorder_pairs)
        bondpairs = set_high_bondorder_pairs(bondpairs, high_bondorder_pairs)

    custom_bondwidth = params.pop('bond_width', 0.2)

    ## C-C triple bond
    #high_bondorder_pairs = {}
    #for i, atom in enumerate(atoms):
    #    if atom.symbol == 'C':
    #        for j in range(i+1,len(atoms)):
    #            if atoms[j].symbol == 'C':
    #                if norm(atoms[i].position-atoms[j].position) < 1.3:
    #                    print(i,j)
    #                    high_bondorder_pairs[(i,j)] = ((0, 0, 0), 3, (0.15, 0.15, 0))

    # set textures 
    textures = []
    for a in atoms:
        textures.append('ase3')

    #colors = []
    #for sym in atoms.get_chemical_symbols():
    #    colors.append(elem_color[sym])
    colors = None
    
    # Common kwargs for eps, png, pov
    kwargs = {
        'rotation'      : custom_rot, # text string with rotation (default='' )
        'radii'         : custom_radii, # float, or a list with one float per atom
        'colors'        : custom_colors,# List: one (r, g, b) tuple per atom
        'show_unit_cell': params['show_box'],   # 0, 1, or 2 to not show, show, and show all of cell
        'celllinewidth' : 0.2,  # Radius of the cylinders representing the cell
        'bondatoms'     : bondpairs, # list of tuples 
        'bondlinewidth' : custom_bondwidth, # linewidth of bond 
        'textures'      : textures, # Length of atoms list of texture names
        }
    
    # Extra kwargs only available for povray (All units in angstrom)
    kwargs.update({
        'run_povray'   : True, # Run povray or just write .pov + .ini files
        'display'      : False,# Display while rendering
        'pause'        : True, # Pause when done rendering (only if display)
        'transparent'  : True,# Transparent background
        'canvas_width' : None, # Width of canvas in pixels
        'canvas_height': 800, # Height of canvas in pixels 
        'camera_dist'  : 50.,  # Distance from camera to front atom
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

    print('Write pov to png.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)

    # positions
    parser.add_argument('params', default='params.json', \
            help='parameter file')

    # options
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-p', '--pos', nargs='?', default='POSCAR.vasp', \
            help='POSCAR File')
    parser.add_argument('-d', '--display', action='store_true', \
            help='if view the structure')

    args = parser.parse_args()

    # setting
    pos_file = args.pos
    name, ext = os.path.splitext(os.path.basename(pos_file))

    with open(args.params, 'r') as reader:
        params = json.load(reader)

    #print(params)

    # render
    atoms = read(pos_file) # 
    #print(atoms)
    if args.display:
        view(atoms)
    else:
        write_pov_png(atoms, name, params)
