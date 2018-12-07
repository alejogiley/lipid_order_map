#!/usr/bin/env python 
##################################################################################################################################################################
#
# This script calculates the order parameters of all
# MARTINI CG lipid C-C bonds with respect to the membrane
# normal. 
#
#	P2 = 0.5 * (3 * <cos^2(theta)> - 1)
#
# where theta is the angle between the bond and the bilayer normal: 
#
# P2 = 1 means perfect alignment with the bilayer normal
# P2 = -0.5 anti-alignment
# P2 = 0 random orientation.
#
# The average vale along all acyl chain bonds
# is stored in a 2D grid, where the bin coordinates are assigned
# according to the (x,y) position of the lipid PO4 bead.  
# Finally, the order-parameter matrix is average over time. 
# ****the membrane normal vector is assumed to be (0,0,1)
#
# hardcode variables:
#
#	stride = 1
#	gridsize = 1
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
import argparse
import numpy as np
import MDAnalysis as mda
import multiprocessing as mp
from scipy import interpolate
from MDAnalysis.analysis.leaflet import LeafletFinder

##################################################################################################################################################################
# use argparse to read in the options from the shell command line
##################################################################################################################################################################

def parse_cmdline(cmdlineArgs):
	"""
	This function initializes the command line parser.
	"""
	
	parser = argparse.ArgumentParser(description='Calculate the order parameter of each acyl-chain bond. \nOutput the time average order parameter in XY plane. \nThe \'trajectory stride\' and \'XY plane gridsize\' are hardcoded',
						usage='use "%(prog)s --help" for more information', 
						formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-p", "--pdbFile", required=True, help='coordinate file: .tpr, .pdb, .gro')
	parser.add_argument("-t", "--trjFile", required=True, help='trajectory file: .dcd, .trr, .xtc')
	parser.add_argument("-l", "--lipidName", help='lipid name: POPE POPG POPC POPS DLPE DLPG DLPC')
	parser.add_argument("-m", "--leafSide", help='monolayer: 1 (up) 2 (down)')
	
	parser.set_defaults(lipidName = "POPC", leafSide = 1)
	
	args = parser.parse_args(cmdlineArgs)

	if (args.pdbFile == None) or (args.trjFile == None):
	
		parser.print_help()
		exit()

	return args.pdbFile, args.trjFile, args.lipidName, args.leafSide

##################################################################################################################################################################
# order function
##################################################################################################################################################################

def Order(number_of_threads, thread_index):
	"""
	This function calculates the cos(theta)^2 and the number of lipids on each grid bin.
	It saves the results in scc.out and cnt.out files, two for each thread.
	File should be combined 'a posteriori' to generate the order_parameter.map
	"""
	
	print "Starting thread {}".format(thread_index)
	
	# read trajectory with mdanalysis
	
	universe = mda.Universe(pdbFile, trjFile)

	# set bead masses

	#universe.atoms[universe.atoms.types == 'G'].masses = 1
	#universe.atoms[universe.atoms.types == 'D'].masses = 1

	# select lipid residues
	
	lipids = universe.select_atoms("resname {}".format(lipidName))

	# save the box dimensions

	(box_x,box_y) = (universe.dimensions[0],universe.dimensions[1])
	
	# save the bead pair for each bond in an array
	
	PO = ['POPE', 'POPG', 'POPC', 'POPS']
	DL = ['DLPE', 'DLPG', 'DLPC']

	if all(x in PO for x in lipidName.split(" ")): 
		bond_names = "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B"
	elif all(x in DL for x in lipidName.split(" ")): 
		bond_names = "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B"
	else:
		print "Error 'lipidName' NOT found in List : " , PO, DL
		os._exit()
	
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
	stride = 1 # analyze every 'stride' frames
	
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
	np.savetxt('tmp_{}/scc_{}_{}.out'.format(leafSide,lipidName[0:2],thread_index), order_array, fmt='%1.4f')
	np.savetxt('tmp_{}/cnt_{}_{}.out'.format(leafSide,lipidName[0:2],thread_index), count_array, fmt='%1.4f')
	
	print "Finished thread {}".format(thread_index)
	
	return 0

##################################################################################################################################################################
# calculate final map
##################################################################################################################################################################

def final(leafSide):
	"""
        This function calculates the final order parameter map
        """
	
	path, dirs, files = next(os.walk('tmp_{}'.format(leafSide)))
	file_count = len(files)
	
	A = 0
	B = 0

	print "Processing files:"

	for name in files:
		print name
		if "cnt" in name:
			a = np.loadtxt('{}/{}'.format(path,name))
			A += a
		elif "scc" in name:
			b = np.loadtxt('{}/{}'.format(path,name))
			B += b

	np.savetxt('test.dat', B, fmt='%1.4f')

	p1 = np.divide(A, B, where=B>=1)
	p2 = 0.5 * (3.0 * p1 - 1.0)

	np.savetxt('order_final_{}.dat'.format(leafSide), p2, fmt='%1.4f')
	
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

	# process the results
	final(leafSide)

	# print the total time
	print("--- %s seconds ---" % (time.time() - start_time))

	exit()

