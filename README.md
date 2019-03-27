# Second-rank order parameter analysis

This script calculates the order parameters of the lipids with respect to the membrane normal.

<p align="center"><img src='./figures/CodeEqn0.gif' /></p>

where theta is the angle between the acyl-chain bond and the bilayer normal, and the average
of <img src='./figures/CodeEqn1.gif'/> is calculated using all the frames in the input ensemble.

```
P2 =  1.0 	perfect alignment with the bilayer normal 
P2 = -0.5 	anti-alignment 
P2 =  0.0 	random orientation.
```

The order parameter of all acyl-chain bonds is stored in a 2D grid, where the bin coordinates 
are assigned according to the x,y position of the lipid's phosphate group.

## Instalation

This script requires:
```
Python 2.7.x
Numpy
Multiprocessing
MDAnalysis
```

Required packages can be installed using [pip](http://www.pip-installer.org/en/latest/index.html), e.g.:
```
pip install MDAnalysis
```

## Citation

Please visit [MDAnalysis](https://www.mdanalysis.org/docs/) website to learn about references.

## Directory tree

```
/msd_map/
├── README.md
├──src
    ├── order_map_parallel.py 	--> python code 
├──files
    ├── system.gro		--> CLC-1 protein in a lipid bilayer
    ├── system.tpr     		--> GRomacs TPR file
    └── traj.trr		--> GRomacs TRR trajectory
└──results 
    ├── order_final_lower.dat	--> average order map for the bottom monolayer
    └── order_final_upper.dat	--> average order map for the top monolayer
```

## Usage

```python
python src/order_map_parallel.py -p files/system.tpr -t files/traj.trr -l POPC -m 1 -a 64 -s 2 
``` 
```
OPTIONS

 -h, --help shows help message and exit

Specify input files and parameters:

 -p	[<.tpr/.pdb/.gro/...>]           (topol.tpr)
          Input coordinate file: .tpr, .pdb, .gro, .psf, ...
 -t	[<.trr/.xtc/.dcd/...>]           (traj.trr)	
	  Input trajectory file: .trr, .xtc, .dcd, .netcdf, .mdcrd, ...
 -l	<string>			 (POPC)	
	  Lipid name: POPE, POPG, POPC, POPS, DLPE, DLPG, DLPC, DPPC
 -m	<int>				 (1)	
	  Monolayer side: 1 (upper) 2 (lower)
 -n	<int>			
	  Number of CPU threads (guess is half of your all your cores)
 -a	<float>				 (8)
	  Area per lipid headgroup in square Angstrom
 -s	<int>				 (1)
	  Sride for the trajectory analysis
```




