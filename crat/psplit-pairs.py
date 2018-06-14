#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 21:06:43 EDT 2018

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
print("[[[[[[[[[[[[[[[[[[[[            PSPLIT-PAIRS                ]]]]]]]]]]]]]]]]]]]]")
print("--------------------------------------------------------------------------------")
print("--------------------------------------------------------------------------------")
print()


# get the environment variable for data directory
ndustdata_dir = os.environ.get('NDUST_DATA')

if not ndustdata_dir:
    ndustdata_dir = '.'

parser = argparse.ArgumentParser()
parser.add_argument("h5prefix", help="h5prefix prefix for h5 and pairs files")
parser.add_argument("-n", "--nchunks", type=int, default=2,
                    help="number of chunks (default: %(default)s)")
parser.add_argument("-d", "--datadir", type=str,
                    default=ndustdata_dir, help="data directory (default: %(default)s)")
# parse arguments
args = parser.parse_args()


h5prefix = args.h5prefix
nchunks = args.nchunks
ndustdata_dir = args.datadir

# change directory
try:
    os.chdir(ndustdata_dir)
except Exception as e:
    print(e)
    sys.exit(-1)

#l = os.listdir()
#print(l)



# split reduced pairs
suffixes = []
try:
    pline = "split --verbose -d -nl/{} --additional-suffix=_rp.dat {}_rp.dat {}-".format(nchunks, h5prefix, h5prefix)
    print(pline)
    p1 = subprocess.Popen(pline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #print('Exit code ', p1.returncode)
    outstd, errstd = p1.communicate()
    suffixes = outstd.decode("utf-8")
    print(suffixes)
    print(errstd.decode("utf-8"))
except Exception as e:
    print(e)
    sys.exit(-1)

# split neutral pair pairs
try:
    pline = "split --verbose -d -nl/{} --additional-suffix=_np.dat {}_np.dat {}-".format(nchunks, h5prefix, h5prefix)
    print(pline)
    p1 = subprocess.Popen(pline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #print('Exit code ', p1.returncode)
    outstd, errstd = p1.communicate()
    print(outstd.decode("utf-8"))
    print(errstd.decode("utf-8"))
except Exception as e:
    print(e)
    sys.exit(-1)

# get the number suffix for each split
numbersplit = [re.sub("_rp.dat'", '', re.sub("creating file '", '', s)) for s in suffixes.splitlines()]

# duplicate h5
for ns in numbersplit:
    try:
        pline = "cp -v {}.h5 {}.h5".format(h5prefix, ns)
        print(pline)
        p1 = subprocess.Popen(pline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #print('Exit code ', p1.returncode)
        outstd, errstd = p1.communicate()
        print(outstd.decode("utf-8"))
        print(errstd.decode("utf-8"))
    except Exception as e:
        print(e)
        sys.exit(-1)

# duplicate particle pairs
for ns in numbersplit:
    try:
        pline = "cp -v {}_pp.dat {}_pp.dat".format(h5prefix, ns)
        print(pline)
        p1 = subprocess.Popen(pline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #print('Exit code ', p1.returncode)
        outstd, errstd = p1.communicate()
        print(outstd.decode("utf-8"))
        print(errstd.decode("utf-8"))
    except Exception as e:
        print(e)
        sys.exit(-1)

# write list of h5 files for crat-merger
filelist = [ns+'.h5' for ns in numbersplit]
with open(h5prefix+".list", mode='wt', encoding='utf-8') as lfile:
    lfile.write('\n'.join(filelist))
