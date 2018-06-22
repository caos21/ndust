#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 10:06:43 EDT 2018

@author: ben

"""
__author__ = "Benjamin Santos"
__copyright__ = "Copyright 2018, Benjamin Santos"
__license__ = "Apache 2.0"
__version__ = "0.1.0"
__email__ = "caos21@gmail.com"
__status__ = "Development"

import os
import sys
import subprocess
import argparse
import re
import h5py
import numpy as np

print()
print("--------------------------------------------------------------------------------")
print("--------------------------------------------------------------------------------")
print("[[[[[[[[[[[[[[[[[[[[            EXTRACT-SLAB                ]]]]]]]]]]]]]]]]]]]]")
print("--------------------------------------------------------------------------------")
print("--------------------------------------------------------------------------------")
print()


# get the environment variable for data directory
ndustdata_dir = os.environ.get('NDUST_DATA')

if not ndustdata_dir:
    ndustdata_dir = '.'

parser = argparse.ArgumentParser()
parser.add_argument("src_h5prefix", help="h5prefix prefix for h5 source file")
parser.add_argument("dest_h5prefix", help="h5prefix prefix for h5 destination file")

parser.add_argument("-d", "--datadir", type=str,
                    default=ndustdata_dir, help="data directory (default: %(default)s)")

parser.add_argument("-l", "--lindex", type=str,
                    help="index of list of indices to extract")

parser.add_argument("-i", "--inputfile", type=str, help="input file for list of radii and charge pairs")


# parse arguments
args = parser.parse_args()

src_h5prefix = args.src_h5prefix
dest_h5prefix = args.dest_h5prefix
lindex = args.lindex
ifile = args.inputfile

# get pairs from file
pairs = []
with open(ndustdata_dir+ifile, mode='r', encoding='utf-8') as lfile:    
    for pair in lfile:
        try:
            diameter, charge = pair.strip('\n').split(',')
            pairs.append({'diameter': diameter, 'charge':charge})
        except Exception as e:
            #print(e)
            pass


class ndust_file(object):
    def __init__(self, src_h5filename, dest_h5filename, pairs=[]):
        '''
        '''
        self.src_h5filename, self.dest_h5filename, self.pairs = src_h5filename, dest_h5filename, pairs
        self.src_h5file = []
        self.dest_h5file = []
        self.groups =  ['Enhancement_factor', 'Barrier_location', 'Barrier_potential',
                        'Coagulation_rate', 'Contact_potential']
        
    def open(self):
        # open 
        self.src_h5file =  h5py.File(self.src_h5filename, 'r')

    def close(self):
        # open
        if self.src_h5file is not None:
            self.src_h5file.close()
        if self.dest_h5file is not None:
            self.dest_h5file.close()
        
    def create(self):
        # open 
        self.dest_h5file =  h5py.File(self.dest_h5filename, 'w')
        
    def clone(self):
        # clone groups        
        self.src_h5file.copy('Charge_sections', self.dest_h5file)
        self.src_h5file.copy('Electrostatic_interaction', self.dest_h5file)
        self.src_h5file.copy('Grid_system', self.dest_h5file)
        self.src_h5file.copy('Volume_sections', self.dest_h5file)

#https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    def find_nearest(self, array, value):
        #array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx, array[idx]

    def clone_pairs(self):
        if self.pairs is not None:
            ## group volumes
            self.gvols = self.src_h5file.get("Volume_sections")
            ## pivots in diameters
            self.dpivots = np.array(self.gvols.get("Diameters"))*1E9

            ## group charges
            self.gchgs = self.src_h5file.get("Charge_sections")
            ## charges
            self.qpivots = np.array(self.gchgs.get("Charges"))

            for group in self.groups:
                ## group Enhancement factor
                self.src_gef = self.src_h5file.get(group)
                if self.src_gef is None:
                    continue

                ## group Enhancement factor
                self.dest_gef = self.dest_h5file.get(group)
                
                ## create Enhancement_factor group
                if self.dest_gef is None:
                    self.dest_h5file.create_group(group)
                    self.dest_gef = self.dest_h5file.get(group)

               
                    for pair in self.pairs:
                       print('Extracting : ', pair)
                       frompair_diameter = float(pair['diameter'])
                       #print(frompair_diameter, 
                       #      self.find_nearest(self.dpivots, frompair_diameter))
                       idx_diameter, diameter = self.find_nearest(self.dpivots, frompair_diameter)

                       frompair_charge = float(pair['charge'])
                       #print(frompair_charge, 
                       #      self.find_nearest(self.qpivots, frompair_charge))
                       idx_charge, charge = self.find_nearest(self.qpivots, frompair_charge)

                       print('Diameter : ', diameter)
                       print('Charge : ', charge)

                       slab_index = (idx_charge+idx_diameter*len(self.qpivots))
                       slab_index_str = str(slab_index)

                       print('Slab : ', slab_index_str)
                       # get slab
                       ef_slab = self.src_gef.get(slab_index_str)

                       # copy dataset
                       self.src_h5file.copy(ef_slab, self.dest_gef, name=slab_index_str)


src_h5filename = ndustdata_dir + src_h5prefix + '.h5'
dest_h5filename = ndustdata_dir + dest_h5prefix + '.h5'

print('Source filename : ', src_h5filename)
print('Destination filename : ', dest_h5filename)

ndf = ndust_file(src_h5filename, dest_h5filename, pairs)

ndf.open()
ndf.create()
ndf.clone()
ndf.clone_pairs()
ndf.close()

