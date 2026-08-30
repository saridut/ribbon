#!/usr/bin/env python

import numpy as np

##Surface triangular mesh file from Gmsh. Must be in Gmsh ver. 2.2.
fn_msh = 'twisted_ribbon.msh'
#
##Pov-Ray include file name.
fn_pov = 'twisted_ribbon.inc'
obj_name = 'TwistedRibbon'
header = 'twisted_ribbon'

##Surface triangular mesh file from Gmsh. Must be in Gmsh ver. 2.2.
#fn_msh = 'coiled_ribbon.msh'
#
##Pov-Ray include file name.
#fn_pov = 'coiled_ribbon.inc'
#obj_name = 'CoiledRibbon'
#header = 'coiled_ribbon'

centered = True #{T, F} Whether to center the object

with open(fn_msh, 'r') as fh:
    while True:
        line = fh.readline().rstrip(' \n')
        if not line:
            break
        if line.startswith('$Nodes'):
            nnodes = int(fh.readline().rstrip(' \n').split()[0])
            nodes = np.zeros((nnodes,3))
            for inode in range(nnodes):
                words = fh.readline().rstrip(' \n').split()
                nodes[inode,0] = float(words[1])
                nodes[inode,1] = float(words[2])
                nodes[inode,2] = float(words[3])
        if line.startswith('$Elements'):
            nelems = int(fh.readline().rstrip(' \n').split()[0])
            elems = np.zeros((nelems,3), dtype=np.int32)
            for ielem in range(nelems):
                words = fh.readline().rstrip(' \n').split()
                elems[ielem,0] = int(words[5]) - 1
                elems[ielem,1] = int(words[6]) - 1
                elems[ielem,2] = int(words[7]) - 1

#Convert to inches from millimeters
#nodes /= 25.4

#Centering 
if centered:
    com = np.mean(nodes,0)
    for inode in range(nnodes):
        nodes[inode,:] -= com

#Calculate normals
face_normals = np.zeros((nelems,3))
node_inc_elems = np.zeros((nnodes,), dtype=np.int32) #Number of incident faces
for ielem in range(nelems):
    elem_nodes = elems[ielem,:]
    #Add up the number of incident faces at each node
    for each in elem_nodes:
        node_inc_elems[each] += 1
    u = nodes[elem_nodes[0],:]
    v = nodes[elem_nodes[1],:]
    w = nodes[elem_nodes[2],:]
    uv = v - u; uv /= np.linalg.norm(uv)
    uw = w - u; uw /= np.linalg.norm(uw)
    normal = np.cross(uv, uw); normal /= np.linalg.norm(normal)
    face_normals[ielem,:] = normal

#Calculate node normals
node_normals = np.zeros((nnodes,3))
for ielem in range(nelems):
    elem_nodes = elems[ielem,:]
    for each in elem_nodes:
        node_normals[each,:] += face_normals[ielem,:]

for inode in range(nnodes):
    normal = node_normals[inode,:]
    node_normals[inode,:] /= np.linalg.norm(normal)

with open(fn_pov, 'w') as fh:
    fh.write('//%s\n'%header)
    fh.write('//\n')
    fh.write('#declare %s = \n'%obj_name)
    fh.write('mesh2 {\n')
    fh.write('  vertex_vectors {\n')
    fh.write('      %i,\n'%nnodes)
    for inode in range(nnodes):
        vx = nodes[inode,0]
        vy = nodes[inode,1]
        vz = nodes[inode,2]
        if inode < (nnodes - 1):
            fh.write('      <%0.12g, %0.12g, %0.12g>,\n'%(vx, vy, vz))
        else:
            fh.write('      <%0.12g, %0.12g, %0.12g> \n'%(vx, vy, vz))
    fh.write('  }\n')
    fh.write('  normal_vectors {\n')
    fh.write('      %i,\n'%nnodes)
    for inode in range(nnodes):
        nx = node_normals[inode,0]
        ny = node_normals[inode,1]
        nz = node_normals[inode,2]
        if inode < (nnodes - 1):
            fh.write('      <%0.12g, %0.12g, %0.12g>,\n'%(nx, ny, nz))
        else:
            fh.write('      <%0.12g, %0.12g, %0.12g> \n'%(nx, ny, nz))
    fh.write('  }\n')
    fh.write('  face_indices {\n')
    fh.write('      %i,\n'%nelems)
    for ielem in range(nelems):
        vi = elems[ielem,0]
        vj = elems[ielem,1]
        vk = elems[ielem,2]
        if ielem < (nelems - 1):
            fh.write('      <%i, %i, %i>,\n'%(vi, vj, vk))
        else:
            fh.write('      <%i, %i, %i> \n'%(vi, vj, vk))
    fh.write('  }\n')
    fh.write('  inside_vector <1, 1, 1>\n')
    fh.write('}\n')
