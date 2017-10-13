# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:34:40 EDT 2017

@author: ben
"""
__author__ = "Benjamin Santos"
__copyright__ = "Copyright 2017, Benjamin Santos"
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "caos21@gmail.com"
__status__ = "Development"

import sys

import numpy as np

import h5py

from scipy import constants

# utils for h5 file construction
#
# ------ H5Writable class ------
class H5Writable(object):
  def __init__(self):
    """ h5path can be a group/path of h5obj
    """
    self.writable_list = []
    self.h5obj = []
    self.h5path = self.h5obj
    #
  def toHDF5(self):
    """ Iterates in writable_list and call toHDF5 methods
    works for DataSet amd Attrib classes
    """
    for writable in self.writable_list:
      writable.toHDF5(self.h5path)
#
# ------ DataSet class ------
class DataSet(object):
    """
    """
    def __init__(self, dsetname, dsetrange, dsettyype, dsetval):
        self.dsetname = dsetname
        self.dsetrange = dsetrange
        self.dsettype = dsettyype
        self.dsetval = dsetval

    """
    """
    def toHDF5(self, h5obj):
        self.h5dset = h5obj.create_dataset(self.dsetname, self.dsetrange,
                                            dtype=self.dsettype)

        self.h5dset[...] =  self.dsetval
#        return self.h5dset

class Attrib(object):
    """
    """
    def __init__(self, attname, attval):
        self.attname = attname
        self.attval = attval
    #
    """
    """
    def toHDF5(self, h5obj):
        self.h5att = h5obj.attrs[self.attname] =  self.attval
    #
    """ Returns value of attribute
    """
    def __call__(self):
        return self.attval
#
