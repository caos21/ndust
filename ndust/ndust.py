# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 14:52:43 EDT 2017

@author: ben

Genrate ui pyuic5 -x ndust.ui -o ndustgui.py
"""
__author__ = "Benjamin Santos"
__copyright__ = "Copyright 2017, Benjamin Santos"
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "caos21@gmail.com"
__status__ = "Development"


import sys
import os
#from PyQt5 import QtGui, QtCore
#from PyQt5.QtCore import SIGNAL
#from PyQt5.QtCore import *
#from PyQt5.QtGui import 
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from ndustgui import Ui_MainWindow

import h5py

#import grid model
import mgrid as mg

#import plasma model
import mplasma as mp

#import plasma model
import mnano as mn

# import hdf5 library wrappers


class Window(QMainWindow, Ui_MainWindow):
  """ Main window
  """
  def __init__(self, parent = None):
    """
    """
    super(Window, self).__init__(parent)
    self.ui = Ui_MainWindow()
    self.ui.setupUi(self)

  # connections
    # grid tab
    grid_save_btn = self.ui.buttonBox_grid_save.button(QtWidgets.QDialogButtonBox.Save)
    grid_save_btn.clicked.connect(self.grid_save)

    self.radiobutton_mpc = self.ui.radioButton_MPC
    self.radiobutton_mpc.toggled.connect(self.toggle_terms)

    self.radiobutton_mpcvdw = self.ui.radioButton_MPCVdW
    self.radiobutton_mpcvdw.toggled.connect(self.toggle_terms)
    self.radiobutton_mpcvdw.toggled.connect(self.toggle_hamaker)
    
    self.radiobutton_ipa = self.ui.radioButton_IPA
    self.radiobutton_ipa.toggled.connect(self.toggle_terms)
    
    self.checkbox_tunnel = self.ui.checkBox_tunnel
    self.checkbox_tunnel.toggled.connect(self.toggle_eaffinity)

    self.checkbox_wnu = self.ui.checkBox_wnu
    self.checkbox_wnu.toggled.connect(self.toggle_nucleation)

    self.checkbox_wsg = self.ui.checkBox_wsg
    self.checkbox_wsg.toggled.connect(self.toggle_sgrowth)

    self.checkbox_wco = self.ui.checkBox_wco
    self.checkbox_wco.toggled.connect(self.toggle_coagulation)

    self.checkbox_wch = self.ui.checkBox_wch
    self.checkbox_wch.toggled.connect(self.toggle_charging)

    # plasma tab
    plasma_save_btn = self.ui.buttonBox_plasma_save.button(QtWidgets.QDialogButtonBox.Save)
    plasma_save_btn.clicked.connect(self.plasma_save)

    # nano tab
    nano_save_btn = self.ui.buttonBox_nano_save.button(QtWidgets.QDialogButtonBox.Save)
    nano_save_btn.clicked.connect(self.nano_save)

    # show window
    self.show()

    self.prefix = ""
    self.dirname = "../data/"

    self.server = [" guillimin:~/duster/results/", " cottos:~/duster/results/"]

  # grid
    # h5 filename for grid
    self.gridfilename = ""

    # container for GridSystem
    self.gsys = []

    # container for volume sections
    self.vsections = []
    self.update_vsections()

    self.spinbox_peakpos = self.ui.spinBox_peakpos
    self.spinbox_peakpos.valueChanged.connect(self.update_peakradius)

    # container for charge sections
    self.qsections = []

    # container for description
    self.gdescription = []

    # container for electrostatic interaction
    self.einteraction = []
    #
  # plasma
    # h5 filename for plasma
    self.plasmafilename = ""

    # container for parameters
    self.parameters = []

    # container for electrons
    self.electrons = []

    # container for ions
    self.ions = []

    # container for metastables
    self.metastables = []

    # container for description
    self.pdescription = []
    #
  # nano
    # h5 filename for nano
    self.nanofilename = ""

    # container for nanoparticle parameters
    self.nanoparticles = []

    # container for rates
    self.rates = []

    # container for time
    self.time = []

    # container for density
    self.density = []
    self.update_peakradius()

    # container for description
    self.ndescription = []
  #
  # WARNING TODO update charge sections too
  def update_vsections(self):
    nvsections = self.ui.spinBox_vsecs.value()
    rmin = le2float(self.ui.lineEdit_minrad)*1e-9

    base = le2float(self.ui.lineEdit_base)
    power = le2float(self.ui.lineEdit_power)

    linear = False
    special = False
    small = False
    medium = False
    big =  False
    if self.ui.checkBox_linear.isChecked():
      linear = True
    if self.ui.checkBox_special.isChecked():
      special = True
      nvsections = 40
      self.ui.spinBox_vsecs.setValue(nvsections)
    if self.ui.checkBox_small.isChecked():
      small = True
      nvsections = 48
      self.ui.spinBox_vsecs.setValue(nvsections)
    if self.ui.checkBox_medium.isChecked():
      medium = True
      nvsections = 48
      self.ui.spinBox_vsecs.setValue(nvsections)
    if self.ui.checkBox_big.isChecked():
      big = True
      nvsections = 48
      self.ui.spinBox_vsecs.setValue(nvsections)      
    grid_type = {'linear':linear, 'special':special, 'small':small, 'medium':medium, 'big':big}


    minvoliface, mg.ifaces, mg.vols, mg.rads, mg.diams = mg.compute_sections(nvsections, rmin,
                                                                             base, power, grid_type)
    vmin = mg.ifaces[0]
    strvmin = "{:.4e}".format(vmin)
    value2le(self.ui.lineEdit_minvol, strvmin)
    vmax = mg.ifaces[-1]
    strvmax = "{:.4e}".format(vmax)
    value2le(self.ui.lineEdit_maxvol, strvmax)
    rmax = mg.rads[-1]*1e9
    strrmax = "{:.4f}".format(rmax)
    value2le(self.ui.lineEdit_maxrad, strrmax)

  def update_peakradius(self):
    """ update peakradius line edit
    """
    peakpos = self.ui.spinBox_peakpos.value()
    if peakpos > len(mg.rads)-1:
      peakpos = len(mg.rads)-1
      self.ui.spinBox_peakpos.setValue(peakpos)
    radius_peakpos = mg.rads[peakpos]*1e9
    strrad = "{:.4f}".format(radius_peakpos)
    value2le(self.ui.lineEdit_radius, strrad)
    return peakpos

  def grid_save(self):
    """ Save grid files
    """      
    print("\n[ii] saving grid file\n")
    self.gridfilename = str(self.ui.lineEdit_grid_save_name.displayText())
    # WARNING must check for empty file FIXME
    #if self.gridfilename <> "":
        #self.prefix += "-"

    h5f = h5py.File(self.dirname+self.gridfilename, "w")
    #h5f.attrs["gtime"] = 0.0

    #profile_rerr = le2float(self.ui.lineEdit_ptol)
    #h5f.attrs["profile_rerr"] = profile_rerr

    # Nanoparticle temperature
    temperature = le2float(self.ui.lineEdit_temp)

    # Nanoparticle mass density
    nmdensity = le2float(self.ui.lineEdit_nmdens)

    # instantiate the group GridSystem
    self.gsys = mg.GridSystem(h5f, temperature, nmdensity)

    # write group
    self.gsys.toHDF5()

    # volume sections
    self.update_vsections()

    # volume sections
    nvsections = self.ui.spinBox_vsecs.value()
    rmin = le2float(self.ui.lineEdit_minrad)*1e-9

    base = le2float(self.ui.lineEdit_base)
    power = le2float(self.ui.lineEdit_power)
    linear = False
    special = False
    small = False
    medium = False
    big =  False
    if self.ui.checkBox_linear.isChecked():
      linear = True
    if self.ui.checkBox_special.isChecked():
      special = True
    if self.ui.checkBox_small.isChecked():
      small = True
    if self.ui.checkBox_medium.isChecked():
      medium = True
    if self.ui.checkBox_big.isChecked():
      big = True
    grid_type = {'linear':linear, 'special':special, 'small':small, 'medium':medium, 'big':big}

    self.vsections = mg.VSections(h5f, nvsections, rmin, base, power, grid_type)

    self.vsections.toHDF5()

    # charge sections
      ##FIXME 
    if self.ui.checkBox_special.isChecked():
      special = True
    max_positive = self.ui.spinBox_maxpos.value()
    max_negative = self.ui.spinBox_maxneg.value()
    if not special:
      nqsections = max_positive + max_negative + 1
      print('not special ', nqsections)
    else:
      #nqsections = 50
      #print('special ', nqsections)
      #max_positive = 5
      #max_negative = 296
      #self.ui.spinBox_maxneg.setValue(max_negative)
      #self.ui.spinBox_maxpos.setValue(max_positive)
      nqsections = 100
      print('special ', nqsections)
      max_positive = 5
      max_negative = 296
      self.ui.spinBox_maxneg.setValue(max_negative)
      self.ui.spinBox_maxpos.setValue(max_positive)

    value2le(self.ui.lineEdit_qsecs, nqsections)
    self.qsections = []
    print('type', type(nqsections))
    if big:
      print('BIG')
      nqsections = 30
      value2le(self.ui.lineEdit_qsecs, nqsections)
      print('type', type(nqsections), nqsections)
      self.qsections = mg.QSections(h5f, nqsections, max_positive, max_negative, grid_type)
    elif small:
      print('SMALL')
      nqsections = 32
      value2le(self.ui.lineEdit_qsecs, nqsections)
      print('type', type(nqsections), nqsections)
      self.qsections = mg.QSections(h5f, nqsections, max_positive, max_negative, grid_type)
    elif medium:
      print('MEDIUM')
      nqsections = 33
      value2le(self.ui.lineEdit_qsecs, nqsections)
      print('type', type(nqsections), nqsections)
      self.qsections = mg.QSections(h5f, nqsections, max_positive, max_negative, grid_type)      
    else:
      self.qsections = mg.QSections(h5f, nqsections, max_positive, max_negative, grid_type)
        
    self.qsections.toHDF5()

    # description for grid
    description_text = qte2string(self.ui.textEdit_grid_desc)

    self.gdescription = mg.Description(h5f, description_text)
    self.gdescription.toHDF5()

    # electrostatic interaction
    multiplier = le2float(self.ui.lineEdit_mult)
    dconstant = le2float(self.ui.lineEdit_die)
    terms = self.ui.spinBox_terms.value()
    hamaker = le2float(self.ui.lineEdit_hamaker)
    vdwradius = le2float(self.ui.lineEdit_vdwradius)    
    #
    method = 0# MPC
    if self.ui.radioButton_IPA.isChecked():
      method = 1# IPA
    elif self.ui.radioButton_coul.isChecked():
      method = 2# Coulomb
    elif self.ui.radioButton_Hybrid.isChecked():
      method = 3# Hybrid
    elif self.ui.radioButton_MPCVdW.isChecked():
      method = 4# MPC + VdW
    #
    self.einteraction = mg.EInteraction(h5f, multiplier, dconstant, method,
                                        terms, hamaker, vdwradius)
    self.einteraction.toHDF5()
    # close file
    h5f.close()
    #
  def toggle_terms(self):
    """ Toggle spin box terms if MPC method is selected
    """
    if self.radiobutton_mpc.isChecked() or self.radiobutton_ipa.isChecked() or self.radiobutton_mpcvdw.isChecked():
      #print("MPC Checked")
      self.ui.label_terms.setEnabled(1)
      self.ui.spinBox_terms.setEnabled(1)
    else:
      #print("MPC UNChecked")
      self.ui.label_terms.setEnabled(0)
      self.ui.spinBox_terms.setEnabled(0)
      #
  def toggle_hamaker(self):
    """ Toggle hamaker and Van der Waals radius line edit if MPC+VdW method is selected
    """
    if self.radiobutton_mpcvdw.isChecked():
      self.ui.label_hamaker.setEnabled(1)
      self.ui.lineEdit_hamaker.setEnabled(1)
      self.ui.label_vdwradius.setEnabled(1)
      self.ui.lineEdit_vdwradius.setEnabled(1)
    else:
      self.ui.label_hamaker.setEnabled(0)
      self.ui.lineEdit_hamaker.setEnabled(0)
      self.ui.label_vdwradius.setEnabled(0)
      self.ui.lineEdit_vdwradius.setEnabled(0)      
    #
  def toggle_eaffinity(self):
    """ Toggle eaffinity box if tunnel is checked
    """
    if self.checkbox_tunnel.isChecked():
      self.ui.label_tunnel.setEnabled(1)
      self.ui.label_eaffinity.setEnabled(1)
      self.ui.lineEdit_eaffinity.setEnabled(1)
    else:
      self.ui.label_tunnel.setEnabled(0)
      self.ui.label_eaffinity.setEnabled(0)
      self.ui.lineEdit_eaffinity.setEnabled(0)
#
  def toggle_nucleation(self):
    """ Toggle nucleation_rate box if nucleation is checked
    """
    if self.checkbox_wnu.isChecked():
      self.ui.label_nucleation.setEnabled(1)
      self.ui.lineEdit_nucleation_rate.setEnabled(1)
    else:
      self.ui.label_nucleation.setEnabled(0)
      self.ui.lineEdit_nucleation_rate.setEnabled(0)
#
  def toggle_sgrowth(self):
    """ Toggle sgrowth_rate box if sgrowth is checked
    """
    if self.checkbox_wsg.isChecked():
      self.ui.label_sgrowth.setEnabled(1)
      self.ui.lineEdit_sgrowth.setEnabled(1)
    else:
      self.ui.label_sgrowth.setEnabled(0)
      self.ui.lineEdit_sgrowth.setEnabled(0)
#
  def toggle_coagulation(self):
    """ Toggle coagulation label if coagulation is checked
    """
    if self.checkbox_wco.isChecked():
      self.ui.label_wco.setEnabled(1)
    else:
      self.ui.label_wco.setEnabled(0)
#
  def toggle_charging(self):
    """ Toggle charging label if charging is checked
    """
    if self.checkbox_wch.isChecked():
      self.ui.label_wch.setEnabled(1)
    else:
      self.ui.label_wch.setEnabled(0)
#

  def plasma_save(self):
    """ Save plasma file
    """      
    print("\n[ii] saving plasma file\n")
    self.plasmafilename = str(self.ui.lineEdit_plasma_save_name.displayText())
    # WARNING must check for empty file FIXME

    h5f = h5py.File(self.dirname+self.plasmafilename, "w")

  # Parameters
    # Fixed plasma check
    pfixed = self.ui.checkBox_pfixed.isChecked()

    # Length
    length = le2float(self.ui.lineEdit_length)

    # Temperature
    temperature = le2float(self.ui.lineEdit_temp_gas)

    # Pressure
    pressure = le2float(self.ui.lineEdit_pressure)

    # instantiate the group Parameters
    self.parameters = mp.Parameters(h5f, pfixed, length, temperature, pressure)

    # fill line edit with gas density
    strng = "{:.4e}".format(self.parameters.neutral_density())
    value2le(self.ui.lineEdit_ng, strng)

    # write group
    self.parameters.toHDF5()

  # Electrons
    # Electron mean energy
    emean = le2float(self.ui.lineEdit_emean)

    # Electron density
    ne = le2float(self.ui.lineEdit_ne)

    # instantiate the group Electrons
    self.electrons = mp.Electrons(h5f, emean, ne)

    # write group
    self.electrons.toHDF5()

  # Ions
    # Ion temperature
    itemp = le2float(self.ui.lineEdit_itemp)

    # Ion density
    ni = le2float(self.ui.lineEdit_ni)

    imass = le2float(self.ui.lineEdit_imass)

    # instantiate the group Ions
    self.ions = mp.Ions(h5f, itemp, ni, imass)

    # write group
    self.ions.toHDF5()

  # Metastables
    # Metastable density
    nm = le2float(self.ui.lineEdit_nm)

    # instantiate the group Electrons
    self.metastables = mp.Metastables(h5f, nm)

    # write group
    self.metastables.toHDF5()

  # Description
    # description for plasma
    description_text = qte2string(self.ui.textEdit_plasma_desc)

    self.pdescription = mp.Description(h5f, description_text)

    # write description
    self.pdescription.toHDF5()

    # close file
    h5f.close()
#
  def nano_save(self):
    """ Save nano file
    """
    print("\n[ii] saving nano file\n")
    self.nanofilename = str(self.ui.lineEdit_nano_save_name.displayText())
    # WARNING must check for empty file FIXME

    h5f = h5py.File(self.dirname+self.nanofilename, "w")

  # Nanoparticles
    # Tunnel current check
    tunnel = self.ui.checkBox_tunnel.isChecked()

    # Electron affinity
    eaffinity = le2float(self.ui.lineEdit_eaffinity)

    # Accomodation factor
    accfactor = le2float(self.ui.lineEdit_acc)

    # instantiate the group Nanoparticles
    self.nanoparticles = mn.Nanoparticles(h5f, tunnel, eaffinity, accfactor)

    # write group
    self.nanoparticles.toHDF5()

  # Rates
    # with nucleation
    wnu = self.ui.checkBox_wnu.isChecked()

    # nucleation rate
    nucleation_rate = le2float(self.ui.lineEdit_nucleation_rate)

    # with surface growth
    wsg = self.ui.checkBox_wsg.isChecked()

    # surface growth rate
    sgrowth_rate = le2float(self.ui.lineEdit_sgrowth)

    # with coagulation
    wco = self.ui.checkBox_wco.isChecked()

    # with charging
    wch = self.ui.checkBox_wch.isChecked()

    # instantiate the group Nanoparticles
    self.rates = mn.Rates(h5f, wnu, nucleation_rate,
                          wsg, sgrowth_rate,
                          wco, wch)

    # write group
    self.rates.toHDF5()

  # Time
    # nanoparticle growth delta time
    ndeltat = le2float(self.ui.lineEdit_ndeltat)

    # nanoparticle charging delta time
    qdeltat = le2float(self.ui.lineEdit_qdeltat)

    # stop time
    tstop = le2float(self.ui.lineEdit_tstop)

    # instantiate the group Time
    self.time = mn.Time(h5f, ndeltat, qdeltat, tstop)

    # write group
    self.time.toHDF5()

  # Density
  # nanoparticle initial density
    indens = le2float(self.ui.lineEdit_indens)

    # nanoparticle qtol
    qtol = le2float(self.ui.lineEdit_qtol)

    # volume sections
    #
    distribution = 0# delta default
    if self.ui.radioButton_step.isChecked():
      distribution = 1# step
    elif self.ui.radioButton_gaussian.isChecked():
      distribution = 2# Gaussian

    peakpos = self.update_peakradius()
    width = self.ui.spinBox_width.value()

    # instantiate the group Time
    self.density = mn.Density(h5f, indens, qtol, distribution, peakpos, width)

    # write group
    self.density.toHDF5()

  # Description
    # description for nano
    description_text = qte2string(self.ui.textEdit_nano_desc)

    self.ndescription = mn.Description(h5f, description_text)

    # write description
    self.ndescription.toHDF5()

    # close file
    h5f.close()
#
def le2int(leobj):
  """ Qt LineEdit to integer
  """
  return int(leobj.displayText())

def le2float(leobj):
  """ Qt LineEdit to float
  """
  # WARNING TODO QValidator or exceptions for inputs
  #try:
    #number = float(leobj.displayText())
  #except Exception:
    #QtGui.QMessageBox.about(self, 'Error','Input can only be a number')
    #pass  
  return float(leobj.displayText())

def value2le(leobj, value):
  """ Set value to Qt LineEdit
  """
  leobj.setText(str(value))

def qte2string(qteobj):
  """ Qt LineEdit to integer
  """
  return str(qteobj.toPlainText())

# instantiation
app = QApplication(sys.argv)
window = Window(None)

sys.exit(app.exec_())
