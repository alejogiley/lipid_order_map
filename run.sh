#!/bin/bash

python2.7 ./src/order_map_parallel.py -p ./input/clc-e1_popeg.tpr -t ./input/clc-e1_popeg.xtc -l 'POPE POPG' -m 1 -n 4 -a 64 -s 1 -o ./output

