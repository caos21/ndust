# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 18:49:48 EDT 2017

@author: ben
"""
__author__ = "Benjamin Santos"
__copyright__ = "Copyright 2016, Benjamin Santos"
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "caos21@gmail.com"
__status__ = "Development"

import sys

import os

import numpy as np

import h5py

# time stamp
import datetime as dt

import mh5utils as mh5u

from scipy import constants

# ------ Parameters class ------
class Parameters(mh5u.H5Writable):
  """ Represents the system global parameters
  """
  def __init__(self, h5obj, pfixed, length, temperature, pressure):
    """
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Parameters")
    #
    # Fixed plasma check
    self.pfixed = mh5u.Attrib("pfixed", int(pfixed))
    #
    # Electrode gap in meters
    self.length = mh5u.Attrib("length", length)
    #
    # Temperature in Kelvin
    self.temperature = mh5u.Attrib("temperature", temperature)
    #
    # Pressure Pascal
    self.pressure = mh5u.Attrib("pressure", pressure)
    #
    ng = (self.pressure() / (constants.Boltzmann * self.temperature()))
    self.neutral_density = mh5u.Attrib("neutral_density", ng)
    #
    self.writable_list = [self.pfixed,
                          self.length,
                          self.temperature,
                          self.pressure,
                          self.neutral_density]
#
# ------ Electrons class ------
class Electrons(mh5u.H5Writable):
  """ Represents the density for electrons
  """
  def __init__(self, h5obj, emean, ne):
    """ Initial values
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Electrons")
    # Set emean electronic mean energy in eV
    self.emean = mh5u.Attrib("emean", emean)
    #
    # Set electronic density 1/m3
    self.ne = mh5u.Attrib("ne", ne)
    #
    self.writable_list = [self.emean, self.ne]
#
# ------ Ions class ------
class Ions(mh5u.H5Writable):
  """ Represents the density for ions
  """
  def __init__(self, h5obj, itemp, ni, imass):
    """ Initial values
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Ions")
    # Set ion temperature in Kelvin
    self.itemp = mh5u.Attrib("itemp", itemp)
    # Set ion density 1/m3s
    self.ni = mh5u.Attrib("ni", ni)
    # Set ion mass in Kg
    self.imass = mh5u.Attrib("imass", imass)
    #
    self.writable_list = [self.itemp, self.ni, self.imass]
#
# ------ Metastables class ------
class Metastables(mh5u.H5Writable):
  """ Represents the density for Metastables
  """
  def __init__(self, h5obj, nm):
    """ Initial values
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Metastables")
    # Set metastable density 1/m3s
    self.nm = mh5u.Attrib("nm", nm)
    #
    self.writable_list = [self.nm]
#
# ------ Description class ------
class Description(mh5u.H5Writable):
  """ Represents a text description for the plasma
  """
  def __init__(self, h5obj, description_text=""):
    """ Initial values
    """
    self.h5obj = h5obj
    # root path
    self.h5path = h5obj
    # date and time
    datetime = str(dt.datetime.now())
    #
    self.description = mh5u.Attrib("Description", description_text)
    #
    self.timestamp = mh5u.Attrib("Timestamp", datetime)
    #
    self.timestamp = mh5u.Attrib("Timestamp", datetime)
    #
    sysinfo = os.uname()
    sysinfostr = ["Sysname",
                  "Nodename",
                  "Release",
                  "Version",  
                  "Machine"]
    self.sysinfo = [None]*5
    for i, si in enumerate(sysinfostr):
      self.sysinfo[i] = mh5u.Attrib(si, str(os.uname()[i]))
    
    # list of lists
    writlist = [[self.description, self.timestamp], self.sysinfo]
    # flatten list
    self.writable_list = [y for x in writlist for y in x]
#
