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
  def __init__(self, h5obj, wnu, nucleation_rate, wsg, sgrowth_rate, wco, wch, wsih4, sih4ratio, sih4nmol):
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
    # With SiH4
    self.wsih4 = mh5u.Attrib("wsih4", int(wsih4))
    # SiH4 : gas ratio
    self.sih4ratio = mh5u.Attrib("sih4ratio", sih4ratio)
    # Number of SiH4 per nucleated particle
    self.sih4nmol = mh5u.Attrib("sih4nmol", sih4nmol)
    # mass of SiH4
    sih4mass = 1.67e-27*(28.+4.)
    self.sih4mass = mh5u.Attrib("sih4mass", sih4mass)
    #
    self.writable_list = [self.wnu,
                          self.nucleation_rate,
                          self.wsg,
                          self.sgrowth_rate,
                          self.wco,
                          self.wch,
                          self.wsih4,
                          self.sih4ratio,
                          self.sih4nmol,
                          self.sih4mass]
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
# ------ Density class ------
class Density(mh5u.H5Writable):
  """ Represents the nanoparticle density parameters
  """
  def __init__(self, h5obj, indens, qtol, distribution, peakpos, width,
               withchargewidth, chargewidth):
    """ Initial values
    """
    # hdf5 file object
    self.h5obj = h5obj
    # hdf5 path
    self.h5path = self.h5obj.create_group("Density")
    # Set nanoparticle initial density
    self.indens = mh5u.Attrib("indens", indens)
    #
    # Set qtol
    self.qtol = mh5u.Attrib("qtol", qtol)
    #
    # Set distribution
    self.distribution = mh5u.Attrib("distribution", distribution)
    #
    # Set peakpos
    self.peakpos = mh5u.Attrib("peakpos", peakpos)
    #
    # Set width (in terms of section number)
    self.width = mh5u.Attrib("width", width)
    #
    # Set chargewidth
    self.chargewidth = mh5u.Attrib("chargewidth", int(withchargewidth))
    self.chargenegwidth = mh5u.Attrib("chargenegwidth", chargewidth["negative"])
    self.chargeposwidth = mh5u.Attrib("chargeposwidth", chargewidth["positive"])

    self.writable_list = [self.indens,
                          self.qtol,
                          self.distribution,
                          self.peakpos,
                          self.width,
                          self.chargewidth,
                          self.chargenegwidth,
                          self.chargeposwidth]
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
