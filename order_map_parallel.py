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
import shutil
import argparse
import numpy as np
import MDAnalysis as mda
import multiprocessing as mp
from MDAnalysis.analysis.leaflet import LeafletFinder

from itertools import islice

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

	parser.add_argument("-p", dest="pdbFile", required=True, help='coordinate file: .tpr, .pdb, .gro')
	parser.add_argument("-t", dest="trjFile", required=True, help='trajectory file: .dcd, .trr, .xtc')
	parser.add_argument("-l", dest="lipidName", help='lipid name: POPE POPG POPC POPS DLPE DLPG DLPC DPPC')
	parser.add_argument("-m", dest="leafSide", help='monolayer: 1 (up) 2 (down)')
	parser.add_argument("-n", dest="nThread", help='number of threads')
	parser.add_argument("-a", dest="Apl", help='area per lipid')
	parser.add_argument("-s", dest="Stride", help='trajectory analysis stride')
	
	parser.set_defaults(lipidName = "POPC", leafSide = 1, Apl = 8, Stride = 1)
	
	args = parser.parse_args(cmdlineArgs)

	if not os.path.exists(args.pdbFile) or not os.path.exists(args.trjFile):
		
		print("ERROR!!!! No such file '{}' or '{}'".format(args.pdbFile, args.trjFile))	
		exit()

	return args.pdbFile, args.trjFile, args.lipidName, args.leafSide, args.nThread, args.Apl, args.Stride

##################################################################################################################################################################
# check-and-balances
##################################################################################################################################################################

def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

def is_empty(obj, namespace):
	if obj:
		return
	else:
		raise ValueError("Something went wrong with {}".format(namestr(obj,namespace)[0]))

##################################################################################################################################################################
# order function
##################################################################################################################################################################

def Order(thread_index):
	"""
	This function calculates the cos(theta)^2 and the number of lipids on each grid bin.
	It saves the results in scc.out and cnt.out files, two for each thread.
	File should be combined 'a posteriori' to generate the order_parameter.map
	"""
	
	print "Starting thread {}".format(thread_index)

	# read trajectory with mdanalysis

	universo = mda.Universe(pdbFile, trjFile)
	
        ###########################################################################
        # INITIALIZE
        ###########################################################################
	
	# SET VARIABLES
	
	stride = int(Stride) # analyze every 'stride' frames
	gridsize = np.sqrt(float(Apl)) # Angstrom
	P2 = 0.0 # order function
	cos_theta = 0.0 # cosine
	
	# INITIALIZE ARRAYS

	(box_x,box_y) = (universo.dimensions[0],universo.dimensions[1])	
        order_array = np.zeros((int((box_x * 2)/gridsize), int((box_y * 2)/gridsize)))
        count_array = np.zeros((int((box_x * 2)/gridsize), int((box_y * 2)/gridsize)))

	# SET LIST OF BONDS FOR LIPID
	
        PO = ['POPE', 'POPG', 'POPC', 'POPS']
        DL = ['DLPE', 'DLPG', 'DLPC']
	DP = ['DPPC']
	
        if all(x in PO for x in lipidName.split(" ")):
                bond_names = "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B"
        elif all(x in DL for x in lipidName.split(" ")):
                bond_names = "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B"
	elif all(x in DP for x in lipidName.split(" ")):
		bond_names = "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B"
        else:
		print "ERROR! ",lipidName," NOT found in List : " , PO, DL, DP
		raise ValueError("Something went wrong with ",lipidName)
	
        enlace = []
        for bond_name in bond_names.split():
                enlace.append(bond_name.split("-"))	

	number_bonds = len(enlace) # number of bonds

	###########################
	
	# SELECT LIPID RESIDUES
	
	lipids = universo.select_atoms("resname {}".format(lipidName))
	is_empty(lipids,locals())

	# SEPARATE LIPIDS BY LEAFLET
	
	L = LeafletFinder(universo, 'name PO4', pbc=True)
	leaflet0 = L.groups(0)
	leaflet1 = L.groups(1)
	
	leaf = {}
	
	po4 = universo.atoms.select_atoms('name PO4')
	bilayerCentre = po4.center_of_geometry()[2]
	
	if leaflet0.centroid()[2] > bilayerCentre:
	        leaf[1] = leaflet0
	        leaf[2] = leaflet1
	else:
	        leaf[2] = leaflet0
	        leaf[1] = leaflet1
	
	monolayer = leaf[int(leafSide)].select_atoms("resname {}".format(lipidName))
	is_empty(monolayer,locals())
	
	###########################################################################
	# START TRAJ CYCLE
	###########################################################################

	blockLength = int(len(universo.trajectory)/int(nThread))
	f = thread_index * blockLength + 1
	l = (thread_index + 1)* blockLength 	

	for ts in universo.trajectory[f:l:stride]:

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
			
			pho = lipids.select_atoms("resnum {} and name PO4".format(res.resnum))
			pho_x = int(math.floor((pho.positions[0][0] + box_x/2)/gridsize))
			pho_y = int(math.floor((pho.positions[0][1] + box_y/2)/gridsize))
			order_array[pho_x,pho_y] += P2 / number_bonds
			count_array[pho_x,pho_y] += 1
			P2 = 0.0 
			cos_theta = 0.0
				
		# end residue cycle
	
	###########################################################################
	# END TRAJ CYCLE
	###########################################################################
	
	print "Saving results in thread {}".format(thread_index)
	np.savetxt('tmp_{}/scc_{}_{}.out'.format(leafSide,lipidName[0:2],thread_index), order_array, fmt='%1.4f')
	np.savetxt('tmp_{}/cnt_{}_{}.out'.format(leafSide,lipidName[0:2],thread_index), count_array, fmt='%1.4f')
	
	print "Finished thread {}".format(thread_index)
	
	return 0

##################################################################################################################################################################
# initialize directory
##################################################################################################################################################################

def inicio():
	"""
	This routine clean and create 
	a temporary directory
	"""
	directory = 'tmp_{}'.format(leafSide)
	
	if os.path.exists(directory):
	       try:
	               shutil.rmtree(directory)
	       except OSError:
	               print ("Deletion of the directory %s failed" % directory)
	       else:
	               print ("Successfully delete the directory %s" % directory)
	
	if not os.path.exists(directory):
	       try:
	               os.makedirs(directory)
	       except OSError:
	               print ("Creation of the directory %s failed" % directory)
	       else:
	               print ("Successfully create the directory %s" % directory)
	return

##################################################################################################################################################################
# process results
##################################################################################################################################################################

def final():
	"""
        This function calculates the final order parameter map
        """
	path, dirs, files = next(os.walk('tmp_{}'.format(leafSide)))
	file_count = len(files)
	
	if file_count > 0:
		A = 0.0
		B = 0.0
		
		print "Processing files:"
		
		for name in files:
			print name
			if "scc" in name:
				a = np.loadtxt('{}/{}'.format(path,name))
				A += np.float64(a)
			elif "cnt" in name:
				b = np.loadtxt('{}/{}'.format(path,name))
				B += np.float64(b)
		
		with np.errstate(divide='ignore',invalid='ignore'): p1 = np.divide(A, B)

		p2 = 0.5 * (3.0 * p1 - 1.0)
		
		np.savetxt('order_final_{}.dat'.format(leafSide), p2, fmt='%1.4f')
	else:
		print('Directory was empty. Exiting')
		exit()
	
	return

##################################################################################################################################################################	
# initialized each thread
##################################################################################################################################################################

def init_worker():
	"""
	This is necessary to be able to interrupt the
	script with CTRL-C (or scancel for that matter).
	It makes sure that the workers ignore SIGINT so
	that any SIGINT sent goes to the master process
	"""
	signal.signal(signal.SIGINT, signal.SIG_IGN)

##################################################################################################################################################################
# main function
##################################################################################################################################################################

def main():
	
	# number of available cpus
	maxcpu = mp.cpu_count()
	
	# the number of threads requested
	jobs = int(nThread)
	if jobs <= maxcpu: maxcpu = jobs

	results = []
	
	print("Running {} threads on {} CPUs".format(jobs,maxcpu))
	
	pool = mp.Pool(maxcpu, init_worker)
	
	for i in range(jobs):
		results.append(pool.apply_async(Order, args=(i, )))

	for r in results: r.get()
	
	try:
		time.sleep(60)
	
	except KeyboardInterrupt:
		print "Caught KeyboardInterrupt, terminating workers"
		pool.terminate()
		pool.join()
		sys.exit(1)
	
	else:
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
	pdbFile, trjFile, lipidName, leafSide, nThread, Apl, Stride = parse_cmdline(sys.argv[1:])
	
	# create output directory
	inicio()
		
	# start the analysis
	main()

	# process the results
	final()

	# print the total time
	print("--- %s seconds ---" % (time.time() - start_time))

	exit()

