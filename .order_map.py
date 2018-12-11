import math
import time
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

from MDAnalysis.analysis.leaflet import LeafletFinder
from matplotlib.patches import Circle

# supress warnings
import warnings
warnings.filterwarnings("ignore")

#
pdbfile = "/u/alejandro/Firefly/2017/Martini/CLC-ec1/monomer/PE_PG/100_0_protein/Production-truncated-files/Truncated-protein/rep0/frame_protmem.gro"

d = {}
for i in range(8):
        d["trjfile" + str(i)] = "/u/alejandro/Firefly/2017/Martini/CLC-ec1/monomer/PE_PG/100_0_protein/Production-truncated-files/Truncated-protein/rep{}/traj_cat_w_ct.xtc".format(i)

universe = mda.Universe(pdbfile, [d['trjfile' + str(i)] for i in range(8)], format='XTC')

lipid_type = "POPE"

if lipid_type == "POPE" or lipid_type == "POPG":
	bond_names = "GL1-C1B GL2-C1A C1A-D2A D2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B"

if lipid_type == "DLPE" or lipid_type == "DLPG":
	bond_names = "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B"

bond = []
for bond_name in bond_names.split():
    bond.append(bond_name.split("-"))

#lipids = universe.select_atoms("resname {}".format(lipid_type))
lipids = universe.select_atoms("resname POPE POPG")

# identify the lipids belonging to each leaf
L = LeafletFinder(universe, 'name PO4', pbc=True)
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)

# upper and lower leaf
leaf = {}

po4 = universe.atoms.select_atoms('name PO4')
bilayerCentre = po4.center_of_geometry()[2]

if leaflet0.centroid()[2] > bilayerCentre:
    leaf[1] = leaflet0
    leaf[2] = leaflet1
else:
    leaf[2] = leaflet0
    leaf[1] = leaflet1
    
del L 
del po4
del leaflet0
del leaflet1

# initialize variables
P = 0.0
P2 = 0.0
step = 1.0
cos_theta = 0.0

# dimensions of the structure file
#(box_x,box_y) = (universe.dimensions[0]+10,universe.dimensions[1]+10)

box_x = box_y = 300.0
# create 2 pixel coordinate mesh grids with a step of 3.0 A
Q = np.mgrid[0:box_x:step, 0:box_y:step]

# initialize 2D grid 
scc = np.zeros((int(box_x/step), int(box_y/step)))
cnt = np.zeros((int(box_x/step), int(box_y/step)))

start_time = time.time()

# Analyze trajectory from t0 to tn every dt -> universe.trajectory[t0:tn:dt]
for ts in universe.trajectory[0::4]:
    
    # Analyze every lipid_type residue
    for res in lipids.residues:
    #for res in leaf[1].select_atoms("resname {}".format(lipid_type)).residues:
        
        # Analyze every C-C bond pair on each acyl chain
        for i in range(len(bond)):
     
            selection = "resnum {} and name {} {}".format(res.resnum, bond[i][0], bond[i][1])
            group = lipids.select_atoms(selection) 
            data = group.positions
            
            cc = (np.concatenate((data[0::2] - data[1::2]), axis=0))
            cc_r = np.sqrt(np.sum(np.power(cc, 2), axis=-1))
            
            cos_theta = cc[..., 2] / cc_r
            P += cos_theta * cos_theta
            
            del data 
            del cc 
            del cc_r
            
        ## end bond cycle
            
        # position of the PO4 bead for residue res
        pho = res.universe.select_atoms("resnum {} and name PO4".format(res.resnum))
        
        # assign residue res to a bin in a 2D XY_grid
        k = int(math.floor((pho.positions[0][0])/step))
        j = int(math.floor((pho.positions[0][1])/step))
        
        # save average P to the assigned grid bin 
        scc[k,j] += P / len(bond)
        cnt[k,j] += 1.0
        
        # reinitialize variables
        P = 0.0 
        cos_theta = 0.0 
        
    ## end residue cycle
    
## end traj cycle

print("--- %s seconds ---" % (time.time() - start_time))

# time average of the P2 values
#old_settings = np.seterr(all='ignore')
array = np.divide(scc, cnt)

# calculate the order parameter
array = 0.5 * (3. * array - 1.)

# position of the interface helix
hlx = res.universe.select_atoms("resnum 177 and name BB")
        
# assign residue res to a bin in a 2D XY_grid
k = int(math.floor((hlx.positions[0][0])/step))
j = int(math.floor((hlx.positions[0][1])/step))

print "protein reference = {}, {}".format(k, j)

np.savetxt('order.out', array, fmt='%1.4e')  # use exponential notation

