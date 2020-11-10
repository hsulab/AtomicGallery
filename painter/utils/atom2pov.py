#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

import os 

import argparse

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

element_radius = {
    'C': 0.5, 'H':0.4, 'O': 0.5, 'Al': 0.75, 'Si': 0.75, 'Cu': 0.85, 
    'Zn': 1.00, 'Cr': 1.00, 'Cd': 1.00
}

def write_pov_png(atoms, rot, sbox, povname):
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

    ## C-C triple bond
    #high_bondorder_pairs = {}
    #for i, atom in enumerate(atoms):
    #    if atom.symbol == 'C':
    #        for j in range(i+1,len(atoms)):
    #            if atoms[j].symbol == 'C':
    #                if norm(atoms[i].position-atoms[j].position) < 1.3:
    #                    print(i,j)
    #                    high_bondorder_pairs[(i,j)] = ((0, 0, 0), 3, (0.15, 0.15, 0))

    ## radii
    rad = [element_radius[at.symbol] for at in atoms]

    #for i, atom in enumerate(atoms):
    #    if atom.symbol != 'C' and atom.symbol != 'H':
    #        if norm(atom.position-cu_pos) > 3.5:
    #            rad[i] = 0.01

    # bond
    bondpairs = get_bondpairs(atoms, radius=0.9)
    #bondpairs = set_high_bondorder_pairs(bondpairs, high_bondorder_pairs)

    # set textures 
    textures = []
    for a in atoms:
        textures.append('ase3')

    elem_color = {
        'Zn': (0.49,0.50,0.69), 
        'Cr': (0.49,0.50,0.69), 
        'O': (1,0.05,0.05), 
    }
    #colors = []
    #for sym in atoms.get_chemical_symbols():
    #    colors.append(elem_color[sym])
    colors = None
    
    # Common kwargs for eps, png, pov
    kwargs = {
        'rotation'      : rot, # text string with rotation (default='' )
        'radii'         : rad, # float, or a list with one float per atom
        'colors'        : colors,# List: one (r, g, b) tuple per atom
        'show_unit_cell': sbox,   # 0, 1, or 2 to not show, show, and show all of cell
        'celllinewidth' : 0.2,  # Radius of the cylinders representing the cell
        'bondatoms'     : bondpairs, # list of tuples 
        # 'bondlinewidth' : 0.2, # linewidth of bond 
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
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-p', '--pos', nargs='?', default='POSCAR.vasp', \
            help='POSCAR File')
    parser.add_argument('-r', '--rot', nargs=3, default=[0,0,0], \
            type=int, help='rotation of image')
    parser.add_argument('-sb', '--showbox', default=0, \
            type=int, help='showbox')
    parser.add_argument('-vis', '--view', default=False, \
            type=bool, help='if view the structure')

    args = parser.parse_args()

    # setting
    # vaspfile = 'CdS-opt.vasp'
    vaspfile = args.pos
    name, ext = os.path.splitext(os.path.basename(vaspfile))
    rot = '{}x,{}y,{}z'.format(*args.rot) # found using ag: 'view -> rotate'

    # render
    atoms = read(vaspfile) # 
    print(atoms)
    if args.view:
        view(atoms)
    else:
        write_pov_png(atoms, rot, args.showbox, name)
