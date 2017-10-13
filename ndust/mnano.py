# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 23:43:51 EDT 2017

@author: ben
"""
__author__ = "Benjamin Santos"
__copyright__ = "Copyright 2017, Benjamin Santos"
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

# ------ Nanoparticles class ------
class Nanoparticles(mh5u.H5Writable):
  """ Represents the Nanoparticles parameters
  """
  def __init__(self, h5obj, tunnel, eaffinity, accfactor):
    """
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Nanoparticles")
    #
    # Fixed tunnel current check
    self.tunnel = mh5u.Attrib("tunnel", int(tunnel))
    #
    # Electron affinity
    self.eaffinity = mh5u.Attrib("eaffinity", eaffinity)
    #
    # Accomodation factor
    self.accfactor = mh5u.Attrib("accfactor", accfactor)
    #
    self.writable_list = [self.tunnel,
                          self.eaffinity,
                          self.accfactor]
#
# ------ Rates class ------
class Rates(mh5u.H5Writable):
  """ Represents the rates nucleation, surface growth, coagulation
      and charging.
  """
  def __init__(self, h5obj, wnu, nucleation_rate, wsg, sgrowth_rate, wco, wch):
    """ Initial values
    """
    # hdf5 file object
    self.h5obj = h5obj
    #
    # hdf5 path
    self.h5path = self.h5obj.create_group("Rates")
    #
    # With nucleation
    self.wnu = mh5u.Attrib("wnu", int(wnu))
    #if not wnu:
      #nucleation_rate = 0.0
    #
    # Nucleation rate
    self.nucleation_rate = mh5u.Attrib("nucleation_rate", nucleation_rate)
    #
    # With surface growth
    self.wsg = mh5u.Attrib("wsg", int(wsg))
    #if not wsg:
      #sgrowth_rate = 0.0
    #
    # Surface growth rate
    self.sgrowth_rate = mh5u.Attrib("sgrowth_rate", sgrowth_rate)
    #
    # With coagulation
    self.wco = mh5u.Attrib("wco", int(wco))
    #
    # With charging
    self.wch = mh5u.Attrib("wch", int(wch))
    #
    self.writable_list = [self.wnu,
                          self.nucleation_rate,
                          self.wsg,
                          self.sgrowth_rate,
                          self.wco,
                          self.wch]
#
# ------ Time class ------
class Time(mh5u.H5Writable):
  """ Represents the time parameters
  """
  def __init__(self, h5obj, ndeltat, qdeltat, tstop):
    """ Initial values
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Time")
    # Set nanoparticle delta t
    self.ndeltat = mh5u.Attrib("ndeltat", ndeltat)
    #
    # Set charging delta t
    self.qdeltat = mh5u.Attrib("qdeltat", qdeltat)
    #
    # Set time stop
    self.tstop = mh5u.Attrib("tstop", tstop)
    #
    self.writable_list = [self.ndeltat,
                          self.qdeltat,
                          self.tstop]
#
# ------ Constants class ------
class Constants(mh5u.H5Writable):
  """ Represents the constant parameters
  """
  def __init__(self, h5obj, indens, qtol):
    """ Initial values
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Constants")
    # Set nanoparticle initial density
    self.indens = mh5u.Attrib("indens", indens)
    #
    # Set qtol
    self.qtol = mh5u.Attrib("qtol", qtol)
    #
    self.writable_list = [self.indens,
                          self.qtol]
#
# ------ Description class ------
class Description(mh5u.H5Writable):
  """ Represents a text description for the Nanoparticles
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
