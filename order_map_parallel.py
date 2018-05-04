#!/usr/bin/env python 
##################################################################################################################################################################
#
# This script calculates the order parameters of all
# MARTINI CG lipid C-C bonds with respect to the membrane
# normal. 
#
#	P2 = 0.5 * (3 * <cos^2(theta)> - 1)
#
# The average vale along all acyl chain bonds
# is stored in a 2D grid, where the bin coordinates are assigned
# according to the (x,y) position of the lipid PO4 bead.  
# Finally, the order-parameter matrix is average over time. 
# ****the membrane normal vector is assumed to be (0,0,1)
#
# 2018 Alejandro Gil-Ley, NHLBI/NIH
#
# MDAnalysis v0.17.0
#
##################################################################################################################################################################

import os
import sys
import math
import time
import signal
import optparse
import numpy as np
import MDAnalysis as mda
import multiprocessing as mp
from MDAnalysis.analysis.leaflet import LeafletFinder

##################################################################################################################################################################
# use argparse to read in the options from the shell command line
##################################################################################################################################################################

def parse_cmdline(cmdlineArgs):
	"""
	This function initializes the command line parser.
	"""
	
	parser = optparse.OptionParser("Usage: order_map_parallel.py -args [options]")
	
	parser.add_option("-p", "--pdb", action="store", dest="pdbFile")	
	parser.add_option("-t", "--trj", action="store", dest="trjFile")
	parser.add_option("-l", "--lipids", action="store", dest="lipidName")
	parser.add_option("-m", "--leaf", action="store", dest="leafSide")

	parser.set_defaults(lipidName = "POPE POPG", leaflet = 1)
	
	opts, args = parser.parse_args(cmdlineArgs)
	pdbFile   = opts.pdbFile
	trjFile   = opts.trjFile
	lipidName = opts.lipidName
	leafSide  = opts.leafSide
	
	# all filenames are required
	if (pdbFile == None) or (trjFile == None):

		parser.print_help()
		exit()
	
	return pdbFile, trjFile, lipidName, leafSide

##################################################################################################################################################################
# order function
##################################################################################################################################################################

def Order(number_of_threads, thread_index):
	"""
	This function calculates the order parameter and the number of lipids on each grid bin.
	It saves the results in scc.out and cnt.out files, two for each thread.
	File should be combined 'a posteriori' to generate the order_parameter.map
	"""
	
	print "Starting thread {}".format(thread_index)
	
	# read trajectory with mdanalysis
	
	universe = mda.Universe(pdbFile, trjFile)
	
	lipids = universe.select_atoms("resname {}".format(lipidName))
	(box_x,box_y) = (universe.dimensions[0],universe.dimensions[1])
	
	# save the bead pair for each bond in an array
	
	if lipidName == "POPE POPG": bond_names = "GL1-C1B GL2-C1A C1A-D2A D2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B"
	elif lipidName == "DLPE DLPG": bond_names = "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B"
	
	enlace = []
	for bond_name in bond_names.split():
		enlace.append(bond_name.split("-"))
	
	# separate lipids by leaflet
	
	L = LeafletFinder(universe, 'name PO4', pbc=True)
	leaflet0 = L.groups(0)
	leaflet1 = L.groups(1)
	
	leaf = {}
	
	po4 = universe.atoms.select_atoms('name PO4')
	bilayerCentre = po4.center_of_geometry()[2]
	
	if leaflet0.centroid()[2] > bilayerCentre:
	        leaf[1] = leaflet0
	        leaf[2] = leaflet1
	else:
	        leaf[2] = leaflet0
	        leaf[1] = leaflet1
	
	monolayer = leaf[int(leafSide)].select_atoms("resname {}".format(lipidName))
	
	# set variables
	
	gridsize = 1 # Angstrom
	P2 = 0.0 # order function
	cos_theta = 0.0 # cosine
		
	number_bonds = len(enlace) # number of bonds
	
	blockLength = int(len(universe.trajectory)/number_of_threads) # length of the trajectory fragments
	start_frame = thread_index * blockLength # start on frame
	end_frame = blockLength * (thread_index + 1) # end on frame
	stride = 4 # analyze every 'stride' frames
	
	# initialize array
	
	order_array = np.zeros((int((box_x)/gridsize), int((box_y)/gridsize)))
	count_array = np.zeros((int((box_x)/gridsize), int((box_y)/gridsize)))
	
	# start traj cycle
	for ts in universe.trajectory[start_frame:end_frame:stride]:
		
		# start residue cycle
		for res in monolayer.residues:
			
			# start bonds cycle
			for i in range(number_bonds):
				
				selection = "resnum {} and name {} {}".format(res.resnum, enlace[i][0], enlace[i][1])
				group = lipids.select_atoms(selection) 
				data = group.positions
				cc = (np.concatenate((data[0::2] - data[1::2]), axis=0))
				cc_r = np.sqrt(np.sum(np.power(cc, 2), axis=-1))  
				cos_theta = cc[..., 2] / cc_r
				P2 += cos_theta * cos_theta
			
			# end bond cycle
			
			pho = res.universe.select_atoms("resnum {} and name PO4".format(res.resnum))
			pho_x = int(math.floor((pho.positions[0][0])/gridsize))
			pho_y = int(math.floor((pho.positions[0][1])/gridsize))
			
			order_array[pho_x,pho_y] += P2 / number_bonds
			count_array[pho_x,pho_y] += 1
			P2 = 0.0 
			cos_theta = 0.0
		
		# end residue cycle
	
	# end traj cycle
	
	print "Saving results in thread {}".format(thread_index)
	np.savetxt('tmp_{}/scc_{}_{}.out'.format(leafSide,lipidName[0:2],thread_index), order_array, fmt='%1.4e')
	np.savetxt('tmp_{}/cnt_{}_{}.out'.format(leafSide,lipidName[0:2],thread_index), count_array, fmt='%1.4e')
	
	print "Finished thread {}".format(thread_index)
	
	return 0


##################################################################################################################################################################	
# initialized each thread
##################################################################################################################################################################

def init_worker():
	signal.signal(signal.SIGINT, signal.SIG_IGN)

##################################################################################################################################################################
# main function
##################################################################################################################################################################

def main():
	
	m = mp.Manager()
	maxcpu = mp.cpu_count()
	jobs = maxcpu / 2
	output = []
	
	print "Initializng {} workers".format(jobs)
	pool = mp.Pool(jobs, init_worker)
	
	print "Starting {} jobs".format(jobs)
	
	for i in range(jobs):
		pool.apply_async(Order, args=(jobs, i, ))
	
	try:
		time.sleep(60)
	
	except KeyboardInterrupt:
		print "Caught KeyboardInterrupt, terminating workers"
		pool.terminate()
		pool.join()
	
	else:
		print "Quitting normally"
		pool.close()
		pool.join()
	
	print 'main process exiting..'


##################################################################################################################################################################
# main routine
##################################################################################################################################################################

if __name__ == "__main__":
	
	# start the time count
	start_time = time.time()
	
	# parse the command line
	pdbFile, trjFile, lipidName, leafSide = parse_cmdline(sys.argv[1:])
	
	# create output directory
	directory = 'tmp_{}'.format(leafSide)
	if not os.path.exists(directory):
		os.makedirs(directory)

	# start the analysis
	main()
	
	# print the total time
	print("--- %s seconds ---" % (time.time() - start_time))

	exit()

