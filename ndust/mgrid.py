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

import os

import numpy as np

import h5py

import mh5utils as mh5u

from scipy import constants

# time stamp
import datetime as dt

# 
def cvol(r):
  """ Volume from radius
  """
  return 4.0 * np.pi * np.power(r, 3) / 3.0
#
def crad(v):
  """ Radius from volume
  """
  return np.cbrt(3.0 * v/(4.0*np.pi))
#
# ------ Grid system class ------
class GridSystem(mh5u.H5Writable):
  """ Represents the system global parameters
  """
  def __init__(self, h5obj, temperature, nmdensity):
    """
    """
    # WARNING FIXME TODO input validation
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Grid_system")
    # nanoparticle temperature
    self.temperature = mh5u.Attrib("Temperature", temperature)
    # nanoparticle mass density
    self.nmdensity = mh5u.Attrib("Mass_density", nmdensity)
    #
    self.writable_list = [self.temperature, self.nmdensity]
#
# ------ Volume Sections class ------
class VSections(mh5u.H5Writable):
  """ Represents the volumes of nanoparticles
  """
  def __init__(self, h5obj, nsections, rmin, base, power):
    """ Initial values
    """
    self.h5obj = h5obj
    #
    self.nsections = mh5u.Attrib("NSections", nsections)
    #
    vmin = cvol(rmin)
    # the minimum volume interface that gives vmin = (v0+v1)/2
    minvoliface = 2.0*vmin/(1.0+base) #vmin = (vo+v1)/2 = vo+base*vo/2
    self.miniface = mh5u.Attrib("Min_interface", minvoliface)
    #
    self.base = mh5u.Attrib("Base", base)
    #
    self.power = mh5u.Attrib("Power", power)
    #
    # Create sections
    self.ifaces = minvoliface*np.power(base, np.arange(0, nsections+1))
    self.interfaces = mh5u.DataSet("Interfaces",
                                    (nsections+1, ),
                                    "f",
                                    self.ifaces)
    #
    self.vols = 0.5*(self.ifaces[1:]+self.ifaces[:-1])
    self.volumes = mh5u.DataSet("Volumes",
                                (nsections, ),
                                "f",
                                self.vols)
    #
    self.rads = crad(self.vols)
    self.radii = mh5u.DataSet("Radii",
                              (nsections, ),
                              "f",
                              self.rads)
    #
    self.diameters = mh5u.DataSet("Diameters",
                                  (nsections, ),
                                  "f",
                                  2.0*self.rads)
    self.h5path = self.h5obj.create_group("Volume_sections")
    self.writable_list = [self.nsections,
                          self.miniface,
                          self.base,
                          self.power,
                          self.interfaces,
                          self.volumes,
                          self.radii,
                          self.diameters]
#
# ------ Charge Sections class ------
class QSections(mh5u.H5Writable):
  """ Represents the charges of nanoparticles
  """
  def __init__(self, h5obj, nsections, maxpositive, maxnegative):
    """ Initial values
    """
    self.h5obj = h5obj
    # Create group
    self.h5path = self.h5obj.create_group("Charge_sections")
    #
    self.nsections = mh5u.Attrib("NSections", nsections)
    #
    self.maxpositive = mh5u.Attrib("Max_positive", maxpositive)
    #
    self.maxnegative = mh5u.Attrib("Max_negative", maxnegative)
    # charges
    charges = np.arange(-maxnegative, maxpositive+1)
    self.charges = mh5u.DataSet("Charges",
                                (nsections, ),
                                "f",
                                charges)    
    #
    self.writable_list = [self.nsections,
                          self.maxpositive,
                          self.maxnegative,
                          self.charges]
#
# ------ Description class ------
class Description(mh5u.H5Writable):
  """ Represents a text description for the grid
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
# ------ Electrostatic Interaction class ------
class EInteraction(mh5u.H5Writable):
  """ Represents the parameters for the electrostatic interaction of
      nanoparticles
  """

  def __init__(self, h5obj, multiplier, dconstant, method, terms=25):    
    """ Initial values
    """
    self.h5obj = h5obj
    # Create group
    self.h5path = self.h5obj.create_group("Electrostatic_interaction")
    #
    self.multiplier = mh5u.Attrib("Multiplier", multiplier)
    #
    self.dconstant = mh5u.Attrib("Dielectric_constant", dconstant)
    #
    self.method = mh5u.Attrib("Method", method)
    #
    if method != 0:# method == 0 is Multipolar Coefficients Potential
      terms = 0
    self.terms = mh5u.Attrib("Terms", terms)
    #
    self.writable_list = [self.multiplier,
                          self.dconstant,
                          self.method,
                          self.terms]
#
