# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 09:16:39 2019

@author: 1203087
"""
import time
import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.lib.distances

import numpy as np
from math import e, sqrt, pi

try:
    import matplotlib

    #matplotlib.use('agg')  # no interactive plotting, only save figures
    #import pylab
    from matplotlib import pyplot as plt
    ## This import registers the 3D projection, but is otherwise unused.
    #from mpl_toolkits.mplot3d import Axes3D
    #from mpl_toolkits.mplot3d.art3d import Line3DCollection # noqa: F401 unused import   
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection # noqa: F401 unused import    
    have_matplotlib = True
except ImportError:
    have_matplotlib = False

import csv

def gaussian(m,s,r):
    gauss = 1/(sqrt(2*pi)*s)*e**(-0.5*((r-m)/s)**2)
    return gauss

def gaussian2(r2):
    gauss = e**(-0.5*r2)
    return gauss
#print(gaussian(0,2,4))
#print(gaussian2(0.0),gaussian2(0.0)/sqrt(2*pi))
#input("Press Enter to continue...")

#MDAnalysis.start_logging()
#start = time.time()                #wall-clock time
#start = time.perf_counter()        #highest available resolution to measure a short duration, include time elapsed during sleep
start = time.process_time()        #sum of the system and user CPU time ,does not include time elapsed during sleep
#file='indentation/indentation/dump/aver50.15.data'
file='dump.010000.data'
simulation=MDAnalysis.Universe(file, atom_style='id type x y z')
all_atoms=simulation.atoms
natom=all_atoms.n_atoms
#print(natom)
posits=all_atoms.positions
#print(posits.ndim, posits.shape, posits.dtype)
#print(posits[0,0],posits[0,1],posits[0,2])
#print(posits[250046,0],posits[250046,1],posits[250046,2])
#globleid=all_atoms.ids
#print(globleid.ndim, globleid.shape, globleid.dtype)
#print(globleid[0],globleid[10])
#input("Press Enter to continue...")

#end = time.time()               #wall-clock time
#end = time.perf_counter()      #highest available resolution to measure a short duration, include time elapsed during sleep
end = time.process_time()      #sum of the system and user CPU time ,does not include time elapsed during sleep
print(end - start)

padding = 4.0 * 2.0
rcut = 3.5

xhi, yhi, zhi = simulation.trajectory.ts.dimensions[:3]
#print( xhi, yhi, zhi )
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
#ax.set_aspect('equal')
ax.set(xlabel='x', ylabel='y', zlabel='z')
v = np.array([[0.0, 0.0, 0.0],
              [xhi, 0.0, 0.0],
              [0.0, yhi, 0.0],
              [0.0, 0.0, zhi],
              [xhi, yhi, 0.0],
              [xhi, 0.0, zhi],
              [0.0, yhi, zhi],
              [xhi, yhi, zhi]])
    
#print((v[:, 0], v[:, 1], v[:, 2]))
ax.scatter3D(v[:, 0], v[:, 1], v[:, 2])
verts = [[v[0],v[1],v[5],v[3]],
         [v[0],v[1],v[4],v[2]],
         [v[0],v[2],v[6],v[3]],
         [v[1],v[4],v[7],v[5]],
         [v[2],v[4],v[7],v[6]],
         [v[3],v[5],v[7],v[6]]]

simulation_box=Poly3DCollection(verts, linewidths=1, edgecolors='k', alpha=.15)
simulation_box.set_facecolor('cyan')
ax.add_collection3d(simulation_box)

xcriteria = "prop x > " + str(padding) + " and " + "prop x < " + str(xhi-padding)
ycriteria = "prop y > " + str(padding) + " and " + "prop y < " + str(yhi-padding)
zcriteria = "prop z > " + str(padding) + " and " + "prop z < " + str(zhi-padding)
#print(xcriteria)
selector = xcriteria + " and " + ycriteria + " and " + zcriteria
#print(selector)

v2 = np.array([[padding, padding, padding],
              [xhi-padding, padding, padding],
              [padding, yhi-padding, padding],
              [padding, padding, zhi-padding],
              [xhi-padding, yhi-padding, padding],
              [xhi-padding, padding, zhi-padding],
              [padding, yhi-padding, zhi-padding],
              [xhi-padding, yhi-padding, zhi-padding]])

ax.scatter3D(v2[:, 0], v2[:, 1], v2[:, 2])
verts2 = [[v2[0],v2[1],v2[5],v2[3]],
         [v2[0],v2[1],v2[4],v2[2]],
         [v2[0],v2[2],v2[6],v2[3]],
         [v2[1],v2[4],v2[7],v2[5]],
         [v2[2],v2[4],v2[7],v2[6]],
         [v2[3],v2[5],v2[7],v2[6]]]

shrink_box=Poly3DCollection(verts2, linewidths=1, edgecolors='b', alpha=.2)
shrink_box.set_facecolor('green')
ax.add_collection3d(shrink_box)

atoms_in_shrinkbox=all_atoms.select_atoms(selector)
#shrinkid=atoms_in_shrinkbox.ids
nshrink=atoms_in_shrinkbox.n_atoms
#print(atoms_in_shrinkbox)
#print(shrinkid[0])
#print(nshrink)

#double confirmed with hand-write function and ovito

#conut = 0                        
#for i in range(1,natom):
#    x = all_atoms[i].position[0]
#    if ( x > padding and x < xhi-padding):
#        conut += 1
#print(conut)
grid_length = 1.0
grid_num = round(rcut/grid_length)
grids = 2*grid_num+1

#Check the grid conpoent works correctly.
#center_atom = atoms_in_shrinkbox[46]
#print(center_atom.position[0], center_atom.position[1])
#intensive = np.zeros( [2*grid_num+1, 2*grid_num+1], dtype=np.float32 )
#print(intensive.shape)
#grid_centers = np.zeros( [2*grid_num+1, 2*grid_num+1, 2], dtype=np.float32 )

#for ix in range(-grid_num,grid_num+1):
#    for iy in range(-grid_num,grid_num+1):        
#        grid_centers[ix+grid_num, iy+grid_num, 0]=center_atom.position[0] + ix*grid_length
#        grid_centers[ix+grid_num, iy+grid_num, 1]=center_atom.position[1] + iy*grid_length
#print(grid_centers, grid_centers.shape, grid_centers.dtype)

f = open('features.csv', 'a', newline='')
writer = csv.writer(f)

f2 = open('check.csv', 'a', newline='')
writer2 = csv.writer(f2)

for i in range(nshrink):
    center_atom = atoms_in_shrinkbox[i]
    print ( "center atom ",atoms_in_shrinkbox.ids[i] , "position ", center_atom.position)
    intensive = np.zeros( [grids, grids, grids], dtype=np.float32 )
    #print(intensive.shape)
    grid_centers = np.zeros( [grids, grids, grids, 3], dtype=np.float32 )
    #print(grid_centers.shape)
    for ix in range(-grid_num,grid_num+1):
        ixx=ix+grid_num
        for iy in range(-grid_num,grid_num+1):
            iyy=iy+grid_num
            for iz in range(-grid_num,grid_num+1):                               
                izz=iz+grid_num
                #print(ix,iy,iz)
                grid_centers[ixx, iyy, izz, 0]=center_atom.position[0] + ix*grid_length
                grid_centers[ixx, iyy, izz, 1]=center_atom.position[1] + iy*grid_length
                grid_centers[ixx, iyy, izz, 2]=center_atom.position[2] + iz*grid_length
    #print(grid_centers, grid_centers.shape, grid_centers.dtype)
    
#    xlb = center_atom.position[0] - padding
#    xhb = center_atom.position[0] + padding
#    ylb = center_atom.position[1] - padding
#    yhb = center_atom.position[1] + padding
#    zlb = center_atom.position[2] - padding
#    zhb = center_atom.position[2] + padding
#    #print ( "pov boundary ", xlb, xhb, ylb, yhb, zlb, zhb )
#    xcriteria = "prop x > " + str(xlb) + " and " + "prop x < " + str(xhb)
#    ycriteria = "prop y > " + str(ylb) + " and " + "prop y < " + str(yhb)
#    zcriteria = "prop z > " + str(zlb) + " and " + "prop z < " + str(zhb)
#    selector = xcriteria + " and " + ycriteria + " and " + zcriteria

#    v3 = np.array([[xlb, ylb, zlb],
#              [xhb, ylb, zlb],
#              [xlb, yhb, zlb],
#              [xlb, ylb, zhb],
#              [xhb, yhb, zlb],
#              [xhb, ylb, zhb],
#              [xlb, yhb, zhb],
#              [xhb, yhb, zhb]])
#
#    ax.scatter3D(v3[:, 0], v3[:, 1], v3[:, 2])
#    verts3 = [[v3[0],v3[1],v3[5],v3[3]],
#              [v3[0],v3[1],v3[4],v3[2]],
#              [v3[0],v3[2],v3[6],v3[3]],
#              [v3[1],v3[4],v3[7],v3[5]],
#              [v3[2],v3[4],v3[7],v3[6]],
#              [v3[3],v3[5],v3[7],v3[6]]]
#
#    shrink_box=Poly3DCollection(verts3, linewidths=1, edgecolors='r', alpha=.3)
#    shrink_box.set_facecolor('m')
#    ax.add_collection3d(shrink_box)
#    
#    
#    ax = fig.add_subplot(122, projection='3d')
#    ax.scatter3D((2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0)
#    label = '(%4.2f, %4.2f, %4.2f)' % (center_atom.position[0],center_atom.position[1],center_atom.position[2])
#    ax.text((2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0, label, 'x', fontsize=8)     
#    ax = fig.gca(projection='3d')
#    ax.set_aspect('equal')
#    ax.voxels(grid_centers[:,:,:,0], facecolors=[0, 1, 1, 0.3], edgecolor='k')
#    ax.set(xlabel='x', ylabel='y', zlabel='z')
#    for ix in range(-grid_num,grid_num+1):
#        for iy in range(-grid_num,-grid_num+1):
#            for iz in range(-grid_num,grid_num+1):
#                x, y, z = grid_centers[ix+grid_num, iy+grid_num, iz+grid_num,:]
#                #label = '(%4.2f, %4.2f, %4.2f)' % (x, y, z)
#                label = '(%4.2f)' % (x)
#                ax.text(ix+grid_num, iy+grid_num, iz+grid_num, label, 'x', fontsize=8)
    
    selector = "point " + str(center_atom.position[0]) + " " + str(center_atom.position[1]) + " " + str(center_atom.position[2]) + " "+ str(padding)
    #print(selector)
    
    atoms_in_pov=all_atoms.select_atoms(selector)
    npov=atoms_in_pov.n_atoms
    #print(atoms_in_pov.ids[39], atoms_in_pov[39].posits)
    #print(npov)
    for j in range(npov):
        neigh_atom = atoms_in_pov[j]
        #print ( "neigh atom ",atoms_in_pov.id[j] , "position ", neigh_atom.position)
        for ix in range(-grid_num,grid_num+1):
            ixx=ix+grid_num
            for iy in range(-grid_num,grid_num+1):
                iyy=iy+grid_num
                for iz in range(-grid_num,grid_num+1):                                      
                    izz=iz+grid_num
                    grid_center=grid_centers[ixx, iyy, izz,:]
                    #print("center of ",ix,iy,iz,"is ", grid_center[:])
                    d2=(grid_center[0]-neigh_atom.position[0])**2+(grid_center[1]-neigh_atom.position[1])**2+(grid_center[2]-neigh_atom.position[2])**2
                    #distance=sqrt(d2)
                    #print("distance to ",ix,iy,iz,"is ",distance)
                    #intensive[ix+grid_num, iy+grid_num, iz+grid_num] += gaussian(0.0,1.0,distance)
                    #print("intensive on ",ix,iy,iz,"is ",gaussian(0.0,1.0,distance))
                    intensive[ixx, iyy, izz] += gaussian2(d2)
                    #print("intensive on ",ix,iy,iz,"is ",gaussian2(d2)/sqrt(2*pi))
                    #print("accumulated intensive on ",ix,iy,iz,"is ",intensive[ix+grid_num, iy+grid_num, iz+grid_num])
                #input("Press Enter to continue...")
#    fig2 = plt.figure()                
#    ax = fig2.add_subplot(111, projection='3d')
#    ax.scatter3D((2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0)
#    label = '(%4.2f, %4.2f, %4.2f)' % (center_atom.position[0],center_atom.position[1],center_atom.position[2])
#    ax.text((2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0,(2.0*float(grid_num)+1.0)/2.0, label, 'x', fontsize=8)     
#    ax = fig2.gca(projection='3d')
#    ax.set_aspect('equal')
#    ax.voxels(grid_centers[:,:,:,0], facecolors=[0, 1, 1, 0.2], edgecolor='k')
#    ax.set(xlabel='x', ylabel='y', zlabel='z')
#    for ix in range(-grid_num,grid_num+1):
#        for iy in range(-grid_num,-grid_num+1):
#            for iz in range(-grid_num,grid_num+1):
#                label = '(%4.2f)' % (intensive[ix+grid_num, iy+grid_num, iz+grid_num])
#                ax.text(ix+grid_num, iy+grid_num, iz+grid_num, label, 'x', fontsize=8) 
#    plt.show()
    intensive /= sqrt(2.0*pi)
    #rowdata=np.append(atoms_in_shrinkbox.ids[i], intensive.flatten())
    rowdata=[atoms_in_shrinkbox.ids[i]]
    rowdata.extend(np.around(intensive, decimals=5).flatten().tolist())
    writer.writerow(rowdata)
    #input("Press Enter to continue...")
    #break
f.close()