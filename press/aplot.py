# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 21:10:56 EST 2016

@author: ben
"""
# Using encoding
# -*- coding: utf-8 -*-
__author__ = "Benjamin Santos"
__copyright__ = "Copyright 2017, Benjamin Santos"
__license__ = "Apache v2.0"
__version__ = "0.1.0"
__email__ = "caos21@gmail.com"
__status__ = "Development"


import sys

import os

import h5py

import numpy as np

import matplotlib

import matplotlib.pyplot as plt

import matplotlib.animation as animation

import matplotlib.mlab as mlab

from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter, LogLocator, LogFormatterExponent

from matplotlib import cm

from matplotlib.colors import LogNorm

import seaborn as sns

current_palette = sns.color_palette()

# Set negative contours to be solid instead of dashed:
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

import matplotlib as mpl
mpl.rcParams['font.size'] = 24
mpl.rcParams['lines.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',
       r'\sisetup{detect-all}',
       r'\usepackage{helvet}',
       r'\usepackage{sansmath}',
       r'\sansmath'] 

import xml.etree.ElementTree as ET

from scipy.integrate import simps

def moments_sum(n, x, y):
  """
  """
  print("mu00 = ", np.sum(n))
  print("mu10 = ", np.sum(x*n), " mu01 = ", np.sum(y*n))
  print("mu11 = ", np.sum(x*y*n))
  print("mu20 = ", np.sum(x*x*n), " mu02 = ", np.sum(y*y*n))
  print("mu21 = ", np.sum(x*x*y*n), " mu12 = ", np.sum(y*y*x*n))
  print("mu30 = ", np.sum(x*x*x*n), " mu03 = ", np.sum(y*y*y*n) )

#print
#print "Moments initial distribution - sum"
#moments_sum(cinitial, X, Y)

def moments_simps(n, x, y):
  """
  """
  print("mu00 = ", simps(simps(n, x), y))
  print("mu10 = ", simps(simps(x*n, x), y), " mu01 = ", simps(simps(y*n, y), x))
  print("mu11 = ", simps(simps(x*n, x)*y, y))
  print("mu20 = ", simps(simps(x*x*n, x), y)," mu02 = ", simps(simps(y*y*n, y),x))
  print("mu21 = ", simps(simps(x*x*n, x)*y, y),
        " mu12 = ", simps(simps(y*y*n, y)*x, x))
  print("mu30 = ", simps(simps(x*x*x*n, x), y),
        " mu03 = ", simps(simps(y*y*y*n, y),x))

class xmlcfg():
  def __init__(self, xmlfname = None):
    # define some default values
    self.tree = None
    self.root = None
    self.figure = "figure"
    self.title = "title"
    self.title_fontsize = 12
    self.tick_fontsize = 8
    self.minortick_fontsize = 6
    self.legend_fontsize = 10
    self.legend_markers = 10

    # find all the axis
    self.axes = None
    self.label_fontsize0 = 10
    self.label_fontsize1 = 10

    self.xmlfname = xmlfname
    if xmlfname is None:
      self.xmlfname = "config.xml"
    else:
      self.parse()

  def load(self, xmlfname):
    self.xmlfname = xmlfname
    if xmlfname is None:
      print("[ee] Missing filename.")

  def parse(self):
    # load the XML tree
    self.tree = ET.parse(self.xmlfname)

    # get the root
    self.root = self.tree.getroot()

    # first validation
    if self.root.tag == "plot":
        #print("[ii] Parsing XML file ", self.xmlfname)

        # get values
        self.figure = self.get_value(self.root, "figure", str)
        self.title = self.get_value(self.root, "title", str)
        self.title_fontsize = self.get_value(self.root, "title_fontsize", int)
        self.tick_fontsize = self.get_value(self.root, "tick_fontsize", int)
        self.minortick_fontsize = self.get_value(self.root,
                                                 "minortick_fontsize", int)
        self.legend_fontsize = self.get_value(self.root, "legend_fontsize",
                                              int)
        self.legend_markers = self.get_value(self.root, "legend_markers", int)

        # find all the axis
        self.axes = self.root.findall('axes')
        self.label_fontsize0 = self.get_value(self.axes[0], "label_fontsize",
                                              int)
        self.label_fontsize1 = self.get_value(self.axes[1], "label_fontsize",
                                              int)
        # find all the results
        self.results = self.root.findall('results')

    else:
        print("[ee] XML Bad format")
        sys.exit(-1)

  def get_value(self, mytree, tag, etype):
      """ Get value from XML element

          :param mytree: the XML node
          :param tag: the name tag of the element
          :param etype: the return type value
          :return: the value corresponding to tag
      """
      value_ = mytree.findall(tag)
      if etype == bool:
          if not value_:
              print("[ii] Can\'t find", tag, "in config file.")
              return etype()
  #            return sys.exit(-1)
          else:
  #            print tag, value_[0].text, str2bool(value_[0].text)
              return str2bool(value_[0].text)
      else:
          if not value_:
              print("[ii] Can\'t find", tag, "in config file.")
              return etype()
  #            return sys.exit(-1)
          else:
              return etype(value_[0].text)

  def str2bool(self, strbool):
      """ Converts from string to bool

          :param v: string
          :return: True or False
      """
      return strbool.lower() in ("yes", "true", "t", "1")

class plot():
  def __init__(self, h5nanoprefix = None, h5gridprefix = None, xc = None):

    self.xc = xc
    if xc is None:
      self.xc = xmlcfg()

    self.defpath = r'/home/ben/git/ndust/data/'
    
    self.h5nanoprefix = h5nanoprefix
    #os.chdir(self.defpath)
    self.h5nanoprefix = self.defpath + self.h5nanoprefix

    # Read nano file
    self.nanofname = self.h5nanoprefix + '.h5'
    self.nanofile =  h5py.File(self.nanofname, 'r')

    self.h5gridprefix = h5gridprefix
    self.h5gridprefix = self.defpath + self.h5gridprefix

    # Read grid file
    self.gridfname = self.h5gridprefix + '.h5'
    self.gridfile =  h5py.File(self.gridfname, 'r')


    ## group volumes
    self.gvols = self.gridfile.get("Volume_sections")
    self.vifaces = np.array(self.gvols.get("Interfaces"))

    ## interfaces in diameters nm
    self.vifaces_diam = np.power(6.0*self.vifaces/np.pi, 1.0/3.0)*1E9

    ## WARNING diameter pivots in nanometres
    self.vpivots = np.array(self.gvols.get("Volumes"))

    ## pivots in diameters
    self.dpivots = np.array(self.gvols.get("Diameters"))*1E9

    ## group charges
    self.gchgs = self.gridfile.get("Charge_sections")

    self.qpivots = np.array(self.gchgs.get("Charges"))

    self.width_vpivots = self.vifaces_diam[1:] - self.vifaces_diam[:-1]
    self.width_qpivots = np.ones_like(self.qpivots)

    # make grid
    self.X, self.Y = np.meshgrid(self.dpivots, self.qpivots)


    # density group
    self.gdensity = self.nanofile.get("Density")
    self.result = np.array(self.gdensity.get("density"))
    self.res_gt_0 = self.result > 0
    self.log10res = np.zeros_like(self.result)
    self.log10res[self.res_gt_0] = np.log10(self.result[self.res_gt_0])

    #self.gtotals = self.nanofile.get("totals")
    #self.total_charge = self.gtotals.attrs["total_charge"]

    #self.qneg = self.gtotals.get("total_negative_charge")

    ## read plasma file
    #self.pgdensity = self.plasmafile.get("density")

    #self.electron = self.pgdensity.get("electron")
    #self.edens = np.array(self.electron.get("new_avg_density"))

    #self.ion = self.pgdensity.get("ion")
    #self.idens = np.array(self.ion.get("new_avg_density"))

    #self.meta = self.pgdensity.get("meta")
    #self.mdens = np.array(self.meta.get("new_avg_density"))

    #self.energy = self.pgdensity.get("energy")
    #self.nrgdens = np.array(self.energy.get("new_avg_density"))

    ## read plasma file
    #self.gmisc = self.plasmafile.get("misc")

    #self.norm = np.array(self.gmisc.get("normalized"))

    #read group electric_background
    #self.eb = self.plasmafile.get("electric_background")

    #self.phi = np.array(self.eb.get("phi_tavg"))
    #self.efield = np.array(self.eb.get("efield_tavg"))

    ##get x pos()
    #self.system = self.plasmafile.get("system")
    #self.xpos = np.array(self.system.get("pos"))

    self.nanofile.close()
    self.gridfile.close()

    self.levelsf = None

    self.extents_linear = [0, len(self.dpivots), 0, len(self.qpivots)]

    self.ddens = np.sum(self.result, axis=0)
    self.cdens = np.sum(self.result, axis=1)

  def print_moments(self):
    print("")
    print("Moments results - sum")
    moments_sum(self.result, self.X, self.Y)
    print("")
    print("Max in section")
    print(self.ddens.argmax())
    print(self.dpivots[self.ddens.argmax()])
    print(np.max(self.ddens))
    ##print "Moments initial distribution - simps"
    ##moments_simps(cinitial, vpivots, qpivots)

  def axis(self, axislabel):
    plt.xlabel(axislabel[0], fontsize=self.xc.label_fontsize0)
    plt.ylabel(axislabel[1], fontsize=self.xc.label_fontsize1)

    plt.xticks(fontsize=self.xc.tick_fontsize)
    # WARNING va='bottom' to correct tick position
    # FIXME TODO get all the ticks and relocate them
    plt.yticks(fontsize=self.xc.tick_fontsize, va='bottom')

    plt.xscale('log')
    #plt.yscale('log')

    plt.xlim((self.dpivots[0], self.dpivots[-1]))
    plt.ylim((self.qpivots[0], self.qpivots[-1]))

  def minorticks(self, fig):
    ax1 = fig.gca()
    ax1.xaxis.set_minor_formatter(FormatStrFormatter("%2.1f"))
    ax1.tick_params(axis='both', which='minor',
                    labelsize=self.xc.minortick_fontsize)


  def plot_b3d(self, data, limits=None, title=r"Density", cmap=cm.seismic,
               axislabel=[r'Diameter $(\text{nm)}$', r'$q(\# e)$', r'$\log N$'],
               vmin=None, vmax=None, savename="bars3d.png"):
    fig = plt.figure(title, figsize=(12, 9))
    fig.suptitle(title, size=self.xc.title_fontsize)

    ax = fig.add_subplot(111, projection='3d')
    plt.tick_params(labelsize=self.xc.tick_fontsize)


    # if limits is None:
    xpos = self.X.flatten()
    ypos = self.Y.flatten()
    elements = (len(self.vpivots)) * (len(self.qpivots))
    zpos = np.zeros(elements)
    data[data<1.0]=0.0
    dz = data.flatten()
    # else:
    #   X, Y = np.meshgrid(self.dpivots[limits[0]:limits[1]],
    #                      self.qpivots[limits[2]:limits[3]])
    #   xpos = X.flatten()
    #   ypos = Y.flatten()
    #   print("len x", len(xpos))
    #   print("len y", len(ypos))
    #   elements = (len(xpos)) * (len(ypos))
    #   zpos = np.zeros(elements)
    #   data[data<1.0]=0.0
    #   dz = data[limits[0]:limits[1], limits[2]:limits[3]].flatten()
    # if limits is not None:
    #   ax.set_xlim3d([limits[0], limits[1]])
    #   ax.set_ylim3d([limits[2], limits[3]])


    dx = 0.75*np.ones_like(zpos)
    dy = dx.copy()

    alpha = np.linspace(0.9, 0.4, len(xpos), endpoint=True)

    colors = cmap(np.linspace(0, 1, len(xpos-1)))
    for i in range(len(xpos)):
      ax.bar3d(xpos[i],ypos[i],zpos[i], dx[i], dy[i], dz[i], alpha=alpha[i],
               color=colors[i], linewidth=0, zsort="average")
    if limits is not None:
      ax.autoscale(False, tight=True)
      ax.set_xlim3d([limits[0], limits[1]])
      ax.set_ylim3d([limits[2], limits[3]])

    ax.set_xlabel("Volume", size=self.xc.tick_fontsize)
    ax.set_ylabel("Charge", size=self.xc.tick_fontsize)
    ax.set_zlabel("$N$", size=self.xc.tick_fontsize)
    # ax.ticklabel_format(axis='z', style='sci', scilimits=(0,2.0E-4))
    #ax.view_init(elev=45, azim=45)
    ax.view_init(elev=12, azim=45)
    plt.show()

  def plot_bars3d(self, limits=None):
    self.plot_b3d(self.log10res, limits)

  def pcolor_dens(self, data, msg="Pcolor", cmap=cm.viridis,
                axislabel=[r'Diameter $(\mathrm{nm)}$', r'$q(\# e)$', r'$\log N$'],
                vmin=None, vmax=None, savename="figpc.png"):
    """
    """
    # Create the figure
    fig = plt.figure(msg, figsize=(12, 9))
    fig.suptitle(msg, size=self.xc.title_fontsize)

    if vmin is None:
      vmin=data.min()*0

    if vmax is None:
      vmax=data.max()

    self.axis(axislabel)

    PC = plt.pcolormesh(self.X, self.Y, data, cmap=cmap,
                        vmin=vmin, vmax=vmax)


    plt.grid()

    # Color bar
    CB = plt.colorbar(PC)
    CB.ax.set_ylabel(axislabel[2], fontsize=self.xc.label_fontsize1)
    CB.ax.tick_params(labelsize=self.xc.tick_fontsize)
    plt.savefig(savename, bbox_inches='tight')
    plt.show()

  def plot_pdens(self):
    self.pcolor_dens(self.log10res, u"Density")

  def plot_fcontours(self, data, msg="Filled contours", cmap=cm.viridis,
                    axislabel=[r'Diameter $(\mathrm{nm)}$', r'$q(\# e)$', r'$\log N$'],
                    savename="figc.png"):
    """
    """
    fig = plt.figure(msg, figsize=(12, 9))
    #fig.suptitle(msg, size=self.xc.title_fontsize)

    if self.levelsf is None:
      self.levelsf = np.arange(1, 18, 2)

    self.axis(axislabel)
    self.minorticks(fig)

    CSF = plt.contourf(self.X, self.Y, data, levels=self.levelsf, cmap=cmap, origin='lower',
                      extent=self.extents_linear)

    CS = plt.contour(self.X, self.Y, data, levels=self.levelsf, origin='lower',
                     colors='k', linewidths=1.5, extent=self.extents_linear)

    plt.clabel(CS, self.levelsf, inline=1, fmt='%d', fontsize=self.xc.tick_fontsize,
              colors='k')

    # Color bar
    CB = plt.colorbar(CSF)
    CB.ax.set_ylabel(axislabel[2], fontsize=self.xc.label_fontsize1)
    CB.ax.tick_params(labelsize=self.xc.tick_fontsize)
    CB.add_lines(CS)

    plt.savefig(savename, bbox_inches='tight')
    plt.show()

  def plot_fdens(self):
    self.plot_fcontours(self.log10res, r'Density',
                        savename="rescf.png")

  def axis1d(self, data, pivots, axislabel, log=True, logy=True):
    plt.xlabel(axislabel[0], fontsize=self.xc.label_fontsize0)
    plt.ylabel(axislabel[1], fontsize=self.xc.label_fontsize1)

    plt.xticks(fontsize=self.xc.tick_fontsize)
    plt.yticks(fontsize=self.xc.tick_fontsize)

    if (log):
      plt.xscale('log')

    if (logy):
      plt.yscale('log')

    plt.xlim((pivots[0], pivots[-1]))
    plt.ylim((1, data.max()))

  def plot_bars(self, pivots, data, width, log=True, msg=u"Size distribution",
                axislabel=[r'$Diameter (\text{nm)}$', r'$N(1/m^3)$'],
                savename="dbars.png", color=current_palette[0], figname="ddiam",
                ylim=None, logy=True):
    """
    """
    # Create the figure
    fig = plt.figure(figname, figsize=(12, 9))
    fig.suptitle(msg, size=self.xc.title_fontsize)

    self.axis1d(data, pivots, axislabel, log=log, logy=logy)
    self.minorticks(fig)

    pbars = plt.bar(pivots, data, width, color=color,
                    edgecolor=color, alpha=0.75, linewidth=3, align='center')

  #  pdens = plt.plot(pivots, data, '-', lw=2.5, color=color)
  # interpolation
    #pivnew = []
    #if log is True:
      #pivnew = np.linspace(pivots.min(), pivots.max(), 100*len(pivots))
      #smooth = spline(vpivots, data, pivnew)
      #fint = interp1d(pivots, data, kind='linear')
      #plt.plot(pivnew, fint(pivnew), lw=2.5)
      #fint = interp1d(pivots, data, kind='slinear')
      #plt.plot(pivnew, fint(pivnew), lw=2.5)
      #fint = interp1d(pivots, data, kind='quadratic')
      #plt.plot(pivnew, fint(pivnew), lw=2.5)
      #fint = interp1d(pivots, data, kind='cubic')
      #plt.plot(pivnew, fint(pivnew), lw=2.5)
      #fint = interp1d(pivots, data, kind='nearest')
      #plt.plot(pivnew, fint(pivnew), lw=2.5)
      #fint = interp1d(pivots, data, kind='zero')
      #plt.plot(pivnew, fint(pivnew), lw=2.5)

    if ylim is not None:
      plt.ylim(ylim)
    plt.savefig(savename, bbox_inches='tight')
    plt.show()

  def plot_diams(self, ylim=None, logy=True):
    self.plot_bars(self.dpivots, self.ddens, self.width_vpivots, ylim=ylim, logy=logy)

  def plot_charges(self, ylim=None, logy=True):
    self.plot_bars(self.qpivots, self.cdens, self.width_qpivots, log=False,
                   msg=u"Charge distribution",
                   axislabel=[r'$q(e)$', r'$N(1/m^3)$'], savename="cbars.png",
                   color=current_palette[1], figname="qbar", ylim=ylim, logy=logy)


  def plot_charges_at_volume(self, sect_indexes, zorders, colors, fills,
                             hatches, alphas, lwidths, xlim=(-20, 2), ylim=(1.0, 1e16)):
    # Create the figure
    savename="cbars-d.png"
    figname="qbar"
    ttitle = "Charge distribution"
    fig = plt.figure(figname, figsize=(12, 9))
    fig.suptitle(ttitle, size=self.xc.title_fontsize)
    axislabel=[r'$q(e)$', r'$N(1/m^3)$']
    plt.xlabel(axislabel[0], fontsize=self.xc.label_fontsize0)
    plt.ylabel(axislabel[1], fontsize=self.xc.label_fontsize1)

    plt.xticks(fontsize=self.xc.tick_fontsize)
    plt.yticks(fontsize=self.xc.tick_fontsize)

    plt.yscale('log')
    plt .ylim(ylim)

    plt.xlim(xlim)

    ax1 = fig.gca()
    ax1.xaxis.set_minor_formatter(FormatStrFormatter("%2.1f"))
    ax1.tick_params(axis='both', which='minor', labelsize=self.xc.minortick_fontsize)

    for i, sect_index in enumerate(sect_indexes):
      charge_dens = self.result[:,sect_index]
      #charge_dens = self.result[sect_index,:]
      #print(self.qpivots)
      diam = self.dpivots[sect_index]
      label = u"$d_{np} = " + str(diam)[:5] + "$ nm"
      #pbars = plt.bar(self.qpivots, charge_dens,
                      #self.width_qpivots, color=colors[i], zorder=zorders[i],
                      #label=label, align='center')
      #pbars = plt.bar(self.qpivots, charge_dens,
                      #1, edgecolor=colors[i], color=colors[i], zorder=zorders[i],
                      #label=label, align='center',fill=fills[i],
                      #hatch=hatches[i], alpha=alphas[i], linewidth=lwidths[i])
      pbars = plt.bar(self.qpivots, charge_dens,
                      1, edgecolor=colors[i], color=colors[i], zorder=zorders[i],
                      label=label, align='center',fill=fills[i],
                      hatch=hatches[i], alpha=alphas[i], linewidth=lwidths[i])
    plt.legend(loc='upper left', shadow=True)
    plt.savefig(savename, bbox_inches='tight')
    plt.show()

  def plot_plasma(self):
    FIG = plt.figure('Densities potential and electric field', figsize=(12, 9))
    FIG.suptitle('Profiles', size=20)

    AX1 = FIG.add_subplot(2, 2, 1)
    #AX1.grid(True)  # Enable Grid Lines
    AX1.set_title('Plasma particle densities')
    plt.xlabel('Position (m)')
    plt.ylabel('Density (1/m3)')
    AX1.grid(True)

    #plt.plot(pos, eDensNew)
    SIZEDPLOT, = AX1.plot(self.xpos, self.edens, 'r-', label='electrons', lw=2.5)
    CHARGDPLOT, = AX1.plot(self.xpos, self.idens, 'b-', label='ions', lw=2.5)
    PMETAS, = AX1.plot(self.xpos, self.mdens, 'g-', label='metastables', lw=2.5)
    PNP, = AX1.plot(self.xpos, -self.norm*self.total_charge, 'm-', label='nanoparticules', lw=2.5)
    #PNP, = AX1.plot(self.xpos, self.norm*self.qneg, 'p-', label='nanoparticules', lw=2.5)

    AX2 = FIG.add_subplot(2, 2, 2)
    AX2.set_title('Energy')
    AX2.grid(True)
    plt.xlabel('Position (m)')
    plt.ylabel('Energy (eV)')

    CHARGDPLOT, = AX2.plot(self.xpos, self.nrgdens/self.edens, 'k-', label='energy', lw=2.5)

    # ***** Eenergy
    #AX3 = FIG.add_subplot(2,2,3)
    #AX3.set_title('Electron Energy Density')
    #AX3.grid(True)
    #plt.xlabel('Position (m)')
    #plt.ylabel('Electron Energy Density (eV/m3)')
    #
    #PEnergy, = AX3.plot([], [], label='pEnergy')

    AX3 = FIG.add_subplot(2, 2, 3)
    AX3.set_title('Electric Field')
    AX3.grid(True)
    plt.xlabel('Position (m)')
    plt.ylabel('Electric Field (kV/m)')

    PEFIELD, = AX3.plot(self.xpos, self.efield[:-1]/1000.0, 'k-', label='efield', lw=2.5)


    AX4 = FIG.add_subplot(2, 2, 4)
    AX4.set_title('Electric Potential')
    AX4.grid(True)
    #AX4.set_ylim(0.0,500.0)
    #AX4.set_ylim(-500.0, 500.0)
    plt.xlabel('Position (m)')
    plt.ylabel('Electric Potential (V)')

    PEVOLT, = AX4.plot(self.xpos, self.phi, 'k-', label='potential', lw=2.5)

    plt.tight_layout()

# ======================= Animation =============================

#==

ADensity = []
DPivots = []
DDens = []
DBoxes = []
QPivots = []
QDens = []
QBoxes = []
PlasmaFile = []
NanoFile = []
WidthDPivots = []
WidthQPivots = []
XMeshDiam = []
YMeshChrgs = [] 

class read_results():
  def __init__(self, nanoh5prefix = None, gridh5prefix = None, plasmah5prefix = None, xc = None, defpath = r'/home/ben/git/ndust/data/'):
    self.nanoh5prefix = None
    self.gridh5prefix = None
    self.plasmah5prefix = None

    self.defpath = r'/home/ben/git/ndust/data/'

    if (nanoh5prefix and gridh5prefix and plasmah5prefix) is None:
      print("Error invalid prefix: nanoh5prefix: {}, gridh5prefix: {}, plasmah5prefix: {}".format(nanoh5prefix, gridh5prefix, plasmah5prefix))
      return None
    else:
      self.nanoh5prefix = nanoh5prefix
      self.gridh5prefix = gridh5prefix
      self.plasmah5prefix = plasmah5prefix

    self.xc = xc
    if xc is None:
      self.xc = xmlcfg()

    #os.chdir(self.defpath)

    # Read nano file
    self.nanofname = self.defpath + self.nanoh5prefix + '.h5'
    global NanoFile
    try:
      NanoFile =  h5py.File(self.nanofname, 'r')
    except:
      print("Error in open nano file: {}".format(nanofname))

    # Read plasma file
    self.gridfname = self.defpath + self.gridh5prefix + '.h5'
    global GridFile
    try: 
      GridFile = h5py.File(self.gridfname, 'r')
    except:
      print("Error in open grid file: {}".format(gridfname))
    
    # Read plasma file
    self.plasmafname = self.defpath + self.plasmah5prefix + '.h5'
    global PlasmaFile
    try:
      PlasmaFile = h5py.File(self.plasmafname, 'r')
    except:
      print("Error in open plasma file: {}".format(plasmafname))

    ## group volumes
    self.gvols = GridFile.get("Volume_sections")
    self.vifaces = np.array(self.gvols.get("Interfaces"))

    ## interfaces in diameters in nanometers
    self.vifaces_diam = np.power(6.0*self.vifaces/np.pi, 1.0/3.0)*1E9

    ## WARNING diameter pivots in nanometres
    self.vpivots = np.array(self.gvols.get("Volumes"))

    ## pivots in diameters
    global DPivots
    DPivots = np.array(self.gvols.get("Diameters"))*1E9

    ## group charges
    self.gchgs = GridFile.get("Charge_sections")

    global QPivots
    QPivots = np.array(self.gchgs.get("Charges"))

    global WidthDPivots
    WidthDPivots = self.vifaces_diam[1:] - self.vifaces_diam[:-1]

    global WidthQPivots
    WidthQPivots = np.ones_like(QPivots)

    # make grid
    global XMeshDiam, YMeshChrgs
    XMeshDiam, YMeshChrgs = np.meshgrid(DPivots, QPivots)

    # density group
    self.gdensity = NanoFile.get("Density")

    global ADensity
    ADensity = np.array(self.gdensity.get("density"))
    # self.res_gt_0 = self.result > 0
    # self.log10res = np.zeros_like(self.result)
    # self.log10res[self.res_gt_0] = np.log10(self.result[self.res_gt_0])

    # self.gtotals = self.nanofile.get("totals")
    # self.total_charge = self.gtotals.attrs["total_charge"]

    # #self.qneg = self.gtotals.get("total_negative_charge")

    # # read plasma file
    # self.pgdensity = self.plasmafile.get("density")

    # self.electron = self.pgdensity.get("electron")
    # self.edens = np.array(self.electron.get("new_avg_density"))

    # self.ion = self.pgdensity.get("ion")
    # self.idens = np.array(self.ion.get("new_avg_density"))

    # self.meta = self.pgdensity.get("meta")
    # self.mdens = np.array(self.meta.get("new_avg_density"))

    # self.energy = self.pgdensity.get("energy")
    # self.nrgdens = np.array(self.energy.get("new_avg_density"))

    # # read plasma file
    # self.gmisc = self.plasmafile.get("misc")

    # self.norm = np.array(self.gmisc.get("normalized"))

    # #read group electric_background
    # self.eb = self.plasmafile.get("electric_background")

    # self.phi = np.array(self.eb.get("phi_tavg"))
    # self.efield = np.array(self.eb.get("efield_tavg"))

    # #get x pos()
    # self.system = self.plasmafile.get("system")
    # self.xpos = np.array(self.system.get("pos"))

    # self.nanofile.close()
    # self.plasmafile.close()

    # self.extents_linear = [0, len(self.dpivots), 0, len(self.qpivots)]

    global DDens
    DDens = np.sum(ADensity, axis=0)

    global QDens
    QDens = np.sum(ADensity, axis=1)


def updateplot(i, profilelist):
    """ Loads a frame i from profilelist
    """

    #print(i)
    ddens = profilelist[0]

    qdens = profilelist[1]

    time = profilelist[2]
    # electric_field = profilelist[3]

    # potential = profilelist[4]

    ## Negative ions profile

    # Set the color to profile
    SIZEDPLOT.set_color('b')
    SIZEDPLOT.set_linewidth(2.0)

    for box, dh in zip(DBoxes, ddens[i]):
        box.set_height(dh)

    #SIZEDPLOT.set_data(DPivots, ddens[i])
    #AX1.text(1e14, 30, 'time')# = %.1f' % pendulum.time_elapsed)
    FIG.suptitle(r'slab \#{}, time {:.3}s'.format(i, time[i]))

    # Allow autoscale in all AXes
    AX1.relim()
    AX1.autoscale_view(True, True, True)


    ## Positive ions profile

    # Set the color to profile
    CHARGDPLOT.set_color('r')
    CHARGDPLOT.set_linewidth(2.0)

    #CHARGDPLOT.set_data(QPivots, qdens[i])
    for box, dh in zip(QBoxes, qdens[i]):
        box.set_height(dh)

    #FSN.bar(QPivots, qdens[i])
    # Allow autoscale in all AXes
    AX2.relim()
    AX2.autoscale_view(True, True, True)

#     PEFIELD.set_data(xpos, electric_field[i][1:])

#     # Set the color to profile
#     PEFIELD.set_color('m')
#     PEFIELD.set_linewidth(2.0)

#     # Allow autoscale only in x
#     AX3.relim()
# #    AX4.autoscale_view(True,scalex=True,scaley=None)
#     AX3.autoscale_view(True, True, True)

#     PEVOLT.set_data(xpos, potential[i])

#     # Set the color to profile
#     PEVOLT.set_color('r')
#     PEVOLT.set_linewidth(2.0)

#     # Allow autoscale only in x
#     AX4.relim()
# #    AX4.autoscale_view(True,scalex=True,scaley=None)
#     AX4.autoscale_view(True, True, True)


# --------------------------- Figure setup ----------------------------------
#
FIG = plt.figure('Densities', figsize=(14, 9))
#FIG.suptitle('Profiles', size=20)

AX1 = FIG.add_subplot(1, 2, 1)
#AX1.grid(True)  # Enable Grid Lines
AX1.set_title('Size distribution')
plt.xlabel('Diameters (nm)')
plt.ylabel('Density (1/m3)')
plt.xscale('log')
plt.yscale('log')
plt.ylim([10, 5e16])
#plt.ylim([10, 1e18])
#plt.xlim([0.75, 80])
plt.xlim([0.75, 100])
#plt.xlim([DPivots[0], DPivots[-1]])

DBARS = []

AX1.grid(True)

#plt.plot(pos, eDensNew)
SIZEDPLOT, = AX1.plot([], [], label='size')

AX2 = FIG.add_subplot(1, 2, 2)
AX2.set_title('Charge distribution')
AX2.grid(True)
plt.xlabel(r'Charges ($e$)')
#plt.ylabel(r'Density ($1/m^3$)')
plt.yscale('log')
plt.ylim([10, 5e16])
plt.xlim([-190, 5])
#plt.xlim([-110, 5])
#plt.xlim([-70, 5])
#plt.ylim([QPivots[0], QPivots[-1]])

CHARGDPLOT, = AX2.plot([], [], label='charges')

# ***** Eenergy
#AX3 = FIG.add_subplot(2,2,3)
#AX3.set_title('Electron Energy Density')
#AX3.grid(True)
#plt.xlabel('Position (m)')
#plt.ylabel('Electron Energy Density (eV/m3)')
#
#PEnergy, = AX3.plot([], [], label='pEnergy')

####AX3 = FIG.add_subplot(2, 2, 3)
####AX3.set_title('Electric Field')
####AX3.grid(True)
####plt.xlabel('Position (m)')
####plt.ylabel('Electric Field (V/m)')

####PEFIELD, = AX3.plot([], [], label='pEnergy')


####AX4 = FIG.add_subplot(2, 2, 4)
####AX4.set_title('Electric Potential')
####AX4.grid(True)
#####AX4.set_ylim(0.0,500.0)
#####AX4.set_ylim(-500.0, 500.0)
####plt.xlabel('Position (m)')
####plt.ylabel('Electric Potential (V)')

####PEVOLT, = AX4.plot([], [], label='PEFIELD')

#plt.tight_layout()
