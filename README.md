# MDAnalysis_order_map

This script calculates the order parameters of all
MARTINI CG lipid C-C bonds with respect to the membrane normal. 

P2 = 0.5 * (3 * <cos^2(theta)> - 1)

where theta is the angle between the bond and the bilayer normal: 

P2 = 1 means perfect alignment with the bilayer normal
P2 = -0.5 anti-alignment
P2 = 0 random orientation.

The average vale along all acyl chain bonds
is stored in a 2D grid, where the bin coordinates are assigned
according to the (x,y) position of the lipid PO4 bead.  
Finally, the order-parameter matrix is average over time. 

# Usage

usage: use "order_map_parallel.py --help" for more information

optional arguments:
  -h, --help    show this help message and exit
  -p PDBFILE    coordinate file: .tpr, .pdb, .gro
  -t TRJFILE    trajectory file: .dcd, .trr, .xtc
  -l LIPIDNAME  lipid name: POPE POPG POPC POPS DLPE DLPG DLPC DPPC
  -m LEAFSIDE   monolayer: 1 (up) 2 (down)
  -n NTHREAD    number of threads
  -a APL        area per lipid
  -s STRIDE     trajectory analysis stride
