import os
import sys
import math
import time
import signal
import argparse
import numpy as np
import MDAnalysis as mda
import multiprocessing as mp

from MDAnalysis.analysis.leaflet import LeafletFinder

# supress warnings
import warnings
warnings.filterwarnings("ignore")

##################################################################################################################################################################

def Order(threads, nthread):

	print "Starting thread {}".format(nthread)
	
	# trajectory ###################################
	
	universe = mda.Universe(options.pdb, options.traj)
	lipids = universe.select_atoms("resname {}".format(options.lipids))
	(box_x,box_y) = (universe.dimensions[0],universe.dimensions[1])
	
	# lipids ######################################
	
	if options.lipids == "POPE POPG": bond_names = "GL1-C1B GL2-C1A C1A-D2A D2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B"
	elif options.lipids == "DLPE DLPG": bond_names = "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B"
	
	enlace = []
	for bond_name in bond_names.split():
		enlace.append(bond_name.split("-"))
	
	# leaflet #####################################
	
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

	monolayer = leaf[int(options.leaf)].select_atoms("resname {}".format(options.lipids))
	
	##########################################
	
	step = 1
	P2 = 0.0
	cos_theta = 0.0
	
	nbond = len(enlace)
	size = len(universe.trajectory)
	block = int(size/threads)
	
	begin = nthread * block
	end = block * (nthread + 1)
	skip = 4
	
	box_x = 340
	box_y = 340

	scc = np.zeros((int((box_x)/step), int((box_y)/step)))
	cnt = np.zeros((int((box_x)/step), int((box_y)/step)))
	
	## start traj cycle
	for ts in universe.trajectory[int(begin):int(end):int(skip)]:
		## start residue cycle
		for res in monolayer.residues:
			## start bond cycle
			for i in range(nbond):
				selection = "resnum {} and name {} {}".format(res.resnum, enlace[i][0], enlace[i][1])
				group = lipids.select_atoms(selection) 
				data = group.positions
				cc = (np.concatenate((data[0::2] - data[1::2]), axis=0))
				cc_r = np.sqrt(np.sum(np.power(cc, 2), axis=-1))  
				cos_theta = cc[..., 2] / cc_r
				P2 += cos_theta * cos_theta
			## end bond cycle
			pho = res.universe.select_atoms("resnum {} and name PO4".format(res.resnum))
			k = int(math.floor((pho.positions[0][0])/step))
			j = int(math.floor((pho.positions[0][1])/step))
			scc[k,j] += P2 / nbond
			cnt[k,j] += 1
			P2 = 0.0 
			cos_theta = 0.0
		## end residue cycle
	## end traj cycle
	
	print "Saving results in thread {}".format(nthread)
	np.savetxt('tmp_{}/scc_{}_{}.out'.format(options.leaf,options.lipids[0:2],nthread), scc, fmt='%1.4e')
	np.savetxt('tmp_{}/cnt_{}_{}.out'.format(options.leaf,options.lipids[0:2],nthread), cnt, fmt='%1.4e')
	
	print "Finished thread {}".format(nthread)
	return

	
##################################################################################################################################################################

#def order_wrapped(a, b):
#	"""A thin wrapper to allow proper tracebacks when things go wrong in a thread"""
#	import traceback
#	try:
#		return Order(a, b)
#	except Exception as e:
#		traceback.print_exc()
#		print('')
#		raise e

def init_worker():
	signal.signal(signal.SIGINT, signal.SIG_IGN)


def main():
	
	m = mp.Manager()
	maxcpu = mp.cpu_count()
	jobs = maxcpu/2
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


## start of main routine ###########################################################################################################################################


if __name__ == '__main__':
	
	start_time = time.time()
	
	# use argparse to read in the options from the shell command line
	parser = argparse.ArgumentParser()
	parser.add_argument("--pdb",default="/u/alejandro/Firefly/2017/Martini/CLC-ec1/monomer/PE_PG/traj_50_50_combine_pbc.pdb",help="the name of the coordinate file to read in (e.g. protein.pdb)")
	parser.add_argument("--traj",default="/u/alejandro/Firefly/2017/Martini/CLC-ec1/monomer/PE_PG/traj_50_50_combine_pbc.dcd",help="the name of the trajectory file to read in (e.g. traj.dcd)")
	parser.add_argument("--lipids",default="POPE POPG",help="lipid types with the same acyl chain lenght and type of atoms (e.g. POPE POPG)")
	parser.add_argument("--leaf",default=1,help="monolayer leaflet to be analyzed: 1 for upper and 2 for lower (e.g. 1)")
	options = parser.parse_args()

	directory = 'tmp_{}'.format(options.leaf)
	
	if not os.path.exists(directory):
		os.makedirs(directory)
	
	main()
	
	print("--- %s seconds ---" % (time.time() - start_time))
	
## end of program ################################################################################################################################################
