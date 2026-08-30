#!/usr/bin/env python

#Converts an xyz file to POV-Ray compatible data file.

import numpy as np

##Input xyz file
lab = 'NPL_ML3_105x10_45_caC8'
fn_xyz = lab + '.xyz'
#
##Pov-Ray include file name.
fn_pov = lab + '.povdat'

centered = True #{T, F} Whether to center the object

with open(fn_xyz, 'r') as fh:
    num_atoms = int(fh.readline().rstrip(' \n').split()[0])
    atom_names = []
    atom_pos = np.zeros((num_atoms,3))
    fh.readline() #Throw away this line
    for iatm in range(num_atoms):
        words = fh.readline().rstrip(' \n').split()
        atom_names.append(words[0])
        atom_pos[iatm,0] = float(words[1])
        atom_pos[iatm,1] = float(words[2])
        atom_pos[iatm,2] = float(words[3])

#Centering 
if centered:
    com = np.mean(atom_pos,0)
    for iatm in range(num_atoms):
        atom_pos[iatm,:] -= com

#Write out
with open(fn_pov, 'w') as fh:
    fh.write('%i,\n'%num_atoms)
    for iatm in range(num_atoms):
        nam = atom_names[iatm]
        x = atom_pos[iatm,0]
        y = atom_pos[iatm,1]
        z = atom_pos[iatm,2]
        if iatm < (num_atoms - 1):
            fh.write('"%s", %0.12g, %0.12g, %0.12g,\n'%(nam, x, y, z))
        else:
            fh.write('"%s", %0.12g, %0.12g, %0.12g \n'%(nam, x, y, z))
