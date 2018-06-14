#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:06:43 EDT 2018

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

print()
print("--------------------------------------------------------------------------------")
print("--------------------------------------------------------------------------------")
print("[[[[[[[[[[[[[[[[[[[[          RUN-PSPLIT-PAIRS              ]]]]]]]]]]]]]]]]]]]]")
print("--------------------------------------------------------------------------------")
print("--------------------------------------------------------------------------------")
print()


# get the environment variable for OMP threads
omp_threads = os.environ.get('OMP_NUM_THREADS')

if not omp_threads:
    omp_threads = 14

# get the environment variable for data directory
ndustdata_dir = os.environ.get('NDUST_DATA')

if not ndustdata_dir:
    ndustdata_dir = '.'

# get the environment variable for build directory
ndustbuild_dir = os.environ.get('NDUST_BUILD')

if not ndustbuild_dir:
    ndustbuild_dir = '/home/ben/git/ndust/crat/build/'

parser = argparse.ArgumentParser()
parser.add_argument("h5list", help="list of h5 files")
parser.add_argument("-t", "--threads", type=int,
                    default=omp_threads, help="Number of OpenMP threads (default: %(default)s)")
parser.add_argument("-d", "--datadir", type=str,
                    default=ndustdata_dir, help="data directory (default: %(default)s)")
parser.add_argument("-b", "--builddir", type=str,
                    default=ndustbuild_dir, help="crat build directory (default: %(default)s)")
# parse arguments
args = parser.parse_args()

h5list = args.h5list
ndustdata_dir = args.datadir
ndustbuild_dir = args.builddir
omp_threads = args.threads

with open(ndustdata_dir+h5list, mode='r', encoding='utf-8') as lfile:
    
    for h5file in lfile:
        h5prefix = re.sub(".h5", '', h5file.strip('\n'))
        #print(h5prefix)
        
        try:
            pline = "time {}./crat-reader -t {} -o {} -d{}".format(ndustbuild_dir, omp_threads, h5prefix, ndustdata_dir)
            print(pline)
            p1 = subprocess.Popen(pline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
            outstd, errstd = p1.communicate()
            print(outstd.decode("utf-8"))
            print(errstd.decode("utf-8"))
        except Exception as e:
            print(e)
            sys.exit(-1)

