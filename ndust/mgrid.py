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
  def __init__(self, h5obj, qpivot, temperature, nmdensity):
    """
    """
    # WARNING FIXME TODO input validation
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Grid_system")
    # pivot in charges
    self.qpivot = mh5u.Attrib("Charge_pivot", int(qpivot))
    # nanoparticle temperature
    self.temperature = mh5u.Attrib("Temperature", temperature)
    # nanoparticle mass density
    self.nmdensity = mh5u.Attrib("Mass_density", nmdensity)
    #
    self.writable_list = [self.qpivot, self.temperature, self.nmdensity]
#
#
def compute_sections(nsections, rmin, base, power, grid_type):
  """ Compute sections
      [rmin] = metre
  """
  if grid_type['linear'] and not grid_type['special']:
    rads = np.linspace(rmin, base*nsections*1e-9, nsections)
    rifaces = rads-base*1e-9*0.5
    rifaces = np.append(rifaces, rads[-1]+base*1e-9*0.5)
    diams = 2.0*rads
    vols = cvol(rads)
    ifaces = cvol(rifaces)
    minvoliface = ifaces[0]
    return minvoliface, ifaces, vols, rads, diams

  if grid_type['special']:
    rlow = np.arange(1, 10.5, 0.5)*1e-9*0.5
    rhigh = np.arange(90, 100.5, 0.5)*1e-9*0.5
    rads = np.concatenate((rlow, rhigh))
    rifaces = rads-base*1e-9*0.5
    rifaces = np.append(rifaces, rads[-1]+base*1e-9*0.5)
    diams = 2.0*rads
    vols = cvol(rads)
    ifaces = cvol(rifaces)
    minvoliface = ifaces[0]
    return minvoliface, ifaces, vols, rads, diams

  if grid_type['small']:
    rlow = np.arange(1, 11.5, 0.25)*1e-9*0.5
    rmed = np.array([45, 50, 55])*1e-9*0.5
    rhigh = np.array([90, 95, 100])*1e-9*0.5
    rads = np.concatenate((rlow, rmed, rhigh))    
    rifaces = rads-base*1e-9*0.5
    rifaces = np.append(rifaces, rads[-1]+base*1e-9*0.5)
    diams = 2.0*rads
    vols = cvol(rads)
    ifaces = cvol(rifaces)
    minvoliface = ifaces[0]
    return minvoliface, ifaces, vols, rads, diams

  if grid_type['medium']:
    rlow = np.array([1, 2, 5, 10])*1e-9*0.5
    rmed = np.arange(45, 55.25, 0.25)*1e-9*0.5
    rhigh = np.array([90, 95, 100])*1e-9*0.5
    rads = np.concatenate((rlow, rmed, rhigh))    
    rifaces = rads-base*1e-9*0.5
    rifaces = np.append(rifaces, rads[-1]+base*1e-9*0.5)
    diams = 2.0*rads
    vols = cvol(rads)
    ifaces = cvol(rifaces)
    minvoliface = ifaces[0]
    return minvoliface, ifaces, vols, rads, diams
  
  if grid_type['big']:
    rlow = np.array([1, 2, 5, 10])*1e-9*0.5
    rmed = np.array([45, 50, 55])*1e-9*0.5
    rhigh = np.arange(90, 100.25, 0.25)*1e-9*0.5
    rads = np.concatenate((rlow, rmed, rhigh))
    rifaces = rads-base*1e-9*0.5
    rifaces = np.append(rifaces, rads[-1]+base*1e-9*0.5)
    diams = 2.0*rads
    vols = cvol(rads)
    ifaces = cvol(rifaces)
    minvoliface = ifaces[0]
    return minvoliface, ifaces, vols, rads, diams
  #
  vmin = cvol(rmin)
  print("radii is {}".format(rmin))
  # the minimum volume interface that gives vmin = (v0+v1)/2
  minvoliface = 2.0*vmin/(1.0+base**power) #vmin = (vo+v1)/2 = vo+base*vo/2
  # Create sections
  ifaces = minvoliface*np.power(base, power*np.arange(0, nsections+1))
  #
  vols = 0.5*(ifaces[1:]+ifaces[:-1])
  rads = crad(vols)
  diams = 2.0*rads
  #
  return minvoliface, ifaces, vols, rads, diams
#
# ------ Volume Sections class ------
class VSections(mh5u.H5Writable):
  """ Represents the volumes of nanoparticles
  """
  def __init__(self, h5obj, nsections, rmin, base, power, grid_type):
    """ Initial values
    """
    self.h5obj = h5obj
    #
    self.nsections = mh5u.Attrib("NSections", nsections)
    #
    minvoliface, self.ifaces, self.vols, self.rads, self.diams = compute_sections(nsections,
                                                                                  rmin,base, power,
                                                                                  grid_type)

    self.miniface = mh5u.Attrib("Min_interface", minvoliface)
    #
    self.base = mh5u.Attrib("Base", base)
    #
    self.power = mh5u.Attrib("Power", power)
    #
    # Create sections
    self.interfaces = mh5u.DataSet("Interfaces",
                                    (nsections+1, ),
                                    "f",
                                    self.ifaces)
    #
    self.volumes = mh5u.DataSet("Volumes",
                                (nsections, ),
                                "f",
                                self.vols)
    #
    self.radii = mh5u.DataSet("Radii",
                              (nsections, ),
                              "f",
                              self.rads)
    #
    self.diameters = mh5u.DataSet("Diameters",
                                  (nsections, ),
                                  "f",
                                  self.diams)
    # WARNING FIXME horrible hack

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
#
# ------ Charge Sections class ------
class QSections(mh5u.H5Writable):
  """ Represents the charges of nanoparticles
  """
  def __init__(self, h5obj, nsections, maxpositive, maxnegative, grid_type):
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
    if not grid_type['special']:
      # charges
      charges = np.arange(-maxnegative, maxpositive+1)
    else:
      # charges
      #charges = np.concatenate((np.arange(-maxnegative, -maxnegative+30),
                                #np.arange(-14, 6)))

      charges = np.concatenate((np.arange(-maxnegative, -maxnegative+30),
                np.arange(-maxnegative+32, -14, 5),
                np.arange(-14, 6)))

    if grid_type['small']:
      charges = np.concatenate((np.array([-235, -225, -215]),
                                np.array([-130, -119, -107]),
                                np.arange(-20, 6)))
      #self.maxpositive = 5
      #self.maxnegative = -235
      print('len charges', len(charges))
      #nsecs = len(charges)
      #self.nsections = mh5u.Attrib("NSections", nsecs)

    if grid_type['medium']:
      charges = np.concatenate((np.array([-235, -225, -215]),
                                np.arange(-130, -106),
                                np.array([-20, -10, -2, -1, 0, 1])))
      #self.maxpositive = 1
      #self.maxnegative = -235
      #nsecs = len(charges)
      #self.nsections = mh5u.Attrib("NSections", nsecs)
      print('len charges', len(charges))
      
    if grid_type['big']:
      charges = np.concatenate((np.arange(-235, -214),
                                np.array([-130, -119, -107]),
                                np.array([-20, -10, -2, -1, 0, 1])))
      #self.maxpositive = 1
      #self.maxnegative = -235
      #nsecs = len(charges)
      #self.nsections = mh5u.Attrib("NSections", nsecs)
      print('len charges', len(charges))

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
    if method == 2:# method == 2 is Coulomb
      terms = 0
    self.terms = mh5u.Attrib("Terms", terms)
    #
    self.writable_list = [self.multiplier,
                          self.dconstant,
                          self.method,
                          self.terms]
#
# ------ van der Waals Interaction class ------
class VDWaals(mh5u.H5Writable):
  """ Represents the parameters for the electrostatic interaction of
      nanoparticles
  """

  def __init__(self, h5obj, vdw=0, cutoff=0.0, bf=0):    
    """ Initial values
    """
    self.h5obj = h5obj
    # Create group
    self.h5path = self.h5obj.create_group("vanderWaals_interaction")
    #
    self.vdw = mh5u.Attrib("On", vdw)
    #
    self.cutoff = mh5u.Attrib("Cutoff", cutoff)
    #
    self.bf = mh5u.Attrib("Brute_force", bf)
    #
    self.writable_list = [self.vdw,
                          self.cutoff,
                          self.bf]
#
