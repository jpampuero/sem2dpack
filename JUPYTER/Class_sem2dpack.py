#!/usr/bin/env python3

"""
@Authors :: Elif Oral & Flomin Tchawe

Class for processing SEM2DPACK outputs.
"""


import sys
sys.path.append('./modules')

import scipy
import numpy as np
import time as t
from matplotlib.path import Path
import matplotlib.patches as pt
import matplotlib.pyplot as plt
import matplotlib as mp
import seaborn as sns
from houches_fb import *
import glob
import fnmatch
from scipy.interpolate import griddata as gd
#from matplotlib.mlab import griddata
import matplotlib.animation as anim
import multiprocessing as mp
import os
import wiggle as wig
from scipy.signal import welch
import scipy.signal as sp
from filters import *
import pandas as pd
import warnings
from math import log10,sin
from matplotlib.colors import LogNorm, Normalize
import imageio, datetime, decimal
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors


warnings.filterwarnings("ignore",category=DeprecationWarning)




# Colormap
def make_colormap():
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    c    = mcolors.ColorConverter().to_rgb
    cc = ['#ffffff', '#dfe6ff', '#a5b5da', '#516b8e', '#c5ce9b',
          '#fcab6c', '#f50000']
    cc1  = np.linspace(0,1,6)
    seq = [c(cc[0]), c(cc[1]), cc1[0],
                      c(cc[1]), c(cc[2]), cc1[1],
                      c(cc[2]), c(cc[3]), cc1[2],
                      c(cc[3]), c(cc[4]), cc1[3],
                      c(cc[4]), c(cc[5]), cc1[4],
                      c(cc[5]), c(cc[6]), cc1[5],
                      c(cc[6])]

    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)
###




def set_style(whitegrid=False, scale=0.85):
  # sns.set(font_scale = 1)
  sns.set_style('white')
  if whitegrid: sns.set_style('whitegrid')
  sns.set_context('talk', font_scale=scale)
#


def make_colors():
  # These are the "Tableau 20" colors as RGB.  
  tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
               (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
               (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
               (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
               (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

  # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
  for i in range(len(tableau20)):
      r, g, b = tableau20[i]
      tableau20[i] = (r / 255., g / 255., b / 255.)

  return tableau20
#


def compute_backbone_curve(gref=0.0, mu=0.0, Nspr=10):
   
    x0 = -6e0
    xu = np.log10(0.1e0)
    dx = (xu-x0)/(Nspr-1)

    gamma  = np.zeros(Nspr)
    G  = np.zeros(Nspr)
    R  = np.zeros(Nspr)

    for i in np.arange(Nspr-1)+1:
        gamma[i] = 10.0**(x0+ dx*(i-1))
        G[i] = 1.0/ (1.0+ abs(gamma[i]/ gref))
        R[i] = G[i]* mu* gamma[i]
    return gamma, R, G
###

def compute_backbone_curve_GGmax1(gref=0.0, mu=0.0, Nspr=10):
   
    x0 = -6e0
    xu = np.log10(0.1e0)
    dx = (xu-x0)/(Nspr-1)

    gamma  = np.zeros(Nspr)
    G  = np.zeros(Nspr)
    R  = np.zeros(Nspr)

    for i in np.arange(Nspr-1)+1:
        gamma[i] = 10.0**(x0+ dx*(i-1))
        G[i] = 1.0
        R[i] = G[i]* mu* gamma[i]
    return gamma, R, G
###




def find_tpeak(time, sr, Vbeta):

  # ATTENTION !!!
  # I ignore the mini pulse with peak amplitude less
  # than sr.max()/ xxx
  # but this may not work for other cases !


  multiple_pulse= False
  t_peak =  time [ np.where( sr == sr.max() ) [0] ]
  ii = np.where ( (time < t_peak) & (sr < Vbeta) )[0]
  
  if len(ii) > 0:  
    ii = ii[::-1]       # in descending order
  else: 
    return t_peak, multiple_pulse

  # Check when v = 0
  tzeros=[]; checkpt=ii[0]
  for i in ii:
    if i != checkpt-1:
      tzeros.append(i)
    checkpt = i

  # Check if any pulse exists before the pulse with max. amplitude
  # excepting the numerical noise
  maxims = []; minilim = sr.max()/6.0  #1.e1
  for t in np.arange(len(tzeros)-1):
    tmax = time [ tzeros[t]   ]
    tmin = time [ tzeros[t+1] ]
    maxi = max ( sr[np.where ( (time < tmax) & (time > tmin)) [0]])
    if maxi > minilim : maxims.append(maxi)

  if len(maxims) > 0:
    tpeaks = [ time [ np.where(sr == maxi)[0] ] for maxi in maxims]
    t_peak = min(tpeaks)
    multiple_pulse = True

  return t_peak, multiple_pulse




# Interpolating SEM results
def interp(x,y,z):
    # Set up a regular grid of interpolation points
    Nx, Ny = 200, 350
    xi, yi = np.linspace(500, 1500, Nx), np.linspace(0, 50, Ny)
    zi = griddata(x, y, z, xi, yi, interp='linear')
    return xi, yi, zi



def read_binary (self, filename, ff):
    dynamic = False # to add later maybe
    with open(filename, 'rb') as fid:
        array = np.fromfile(fid, ff) 
        output = np.zeros( (self.npts,self.nsta_extra) )
        if dynamic:
            for i in np.arange(len(array)/self.nsta_extra):
                limit1 = i*self.nsta_extra
                limit2 = (i+1)*self.nsta_extra
                output[i,:] = array[limit1:limit2]
        else:
            for i in np.arange(self.nsta_extra):
                limit1 = i* self.npts
                limit2 = (i+1)* self.npts
                output[:,i] = array[limit1:limit2]
        fid.close()   
    return output
###

def interpg(field,coord, inc=100):
    """
    Interpolates argument field over a meshgrid.
    Meshgrid size depends on argument coord.
    """
    print ('Interpolating...')
    xcoord = coord[:,0]
    zcoord = coord[:,1]
    nbx = len(xcoord)
    nbz = len(zcoord)
    ext = [min(xcoord), max(xcoord), min(zcoord), max(zcoord)]
    x,z=np.meshgrid(np.linspace(ext[0],ext[1],inc),np.linspace(ext[2],ext[3],inc),sparse=True)
    y = gd((xcoord,zcoord),field,(x,z),method='linear')
    y =np.flipud(y)
    return y
#

def create_gif(filenames, duration):
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
    output_file = 'animated_snapshots.gif'
    imageio.mimsave(output_file, images, duration=duration)
    print ('Animation saved as ', output_file)
#


def write_station_file(filename='stats.txt',lims=[0, 1, 0, 1], dx=1.0, dz=1.0):
    f = open (filename, 'w')
    nsta = 0 
    for _z in np.arange(lims[2],lims[3]+dz, dz):
        for _x in np.arange(lims[0],lims[1]+dx, dx):
            f.write('%s \t %s \n' % (str(_x), str(_z)) )
            nsta+=1
    f.close()
    print ('Number of stations (input file): ', nsta)
    return
#



class sem2dpack(object):
  """
    Class to postprocess SEM2DPACK simulation outputs.
    It implements both object methods and static methods
    for easy handling and visulazing the data outputs.
    It has the following instances:
      - directory :: the simulation directory
      - mdict     :: dictionary containing spectral element grid infos
      - dt        :: simulation time step
      - npts      :: Number of points in record, npts * dt gives
               the simulation duration
      - nsta      :: number of reciever stations
      - velocity  :: velocity traces
      - time      :: time vector (0:dt:npts*dt)
      - fmax      :: maximum frequency of simulation
      - tf        :: transfer function in case of sedimentary basins
      - f         :: frequecy vector
      - rcoord    :: reciever stations coordinates
  """

  def __init__(self,directory,freqs=(0.1,10), extra=False, db_precision=False, read_header=True):
    self.directory = directory+ '/'
    self.mdict = {}
    self.dt = 0.0
    self.npts = 0
    self.nsta = 0
    self.velocity = np.array([])
    self.time = np.array([])
    self.fmax = 0.0
    self.tf = np.array([])
    self.f_interp = np.array([])
    self.f = np.array([])
    self.rcoord = np.array([])
    self.extra_coord = np.array([])  
    self.nonlinear_curve = np.array([])
    self.x_interp = np.array([])
    self.vs_int = np.array([])
    self.extra = extra
    self.double_precision = db_precision
    self.Elastic = True
    self.Dynamic = False
    self.Effective = False
    self.fault = {}

    try:
      self.__readSpecgrid()
    except:
      raise Exception('Check directory name \
         Problem reading specgrid files!')
    #

    if read_header:
      try:
        self.__read_header()
      except:
        raise Exception('Check directory name \
           OR no SeisHeader_sem2d.hdr file!')
    #
    print ('*')
    ###
    



#

  def __readSpecgrid(self):
    print ('Reading grid information...')
    # read mesh information
    filename = self.directory + '/grid_sem2d.hdr'
    line = np.genfromtxt(filename, dtype=int)
    nel, npgeo, ngnod, npt, ngll = line[1,:]

    # read spectral element grid coordinates
    fname = self.directory + '/coord_sem2d.tab'
    line = pd.read_csv(fname, delim_whitespace=True, header=None, nrows=1, dtype=int)
    data = pd.read_csv(fname, names=('x','z'), delim_whitespace=True, header=0, nrows=line.values[0][0])
    coord = np.vstack ( (data['x'].values, data['z'].values)  ).T

    #read ibool file
    filename = self.directory + '/ibool_sem2d.dat'
    with open(filename,'rb') as f:
      ibool = np.fromfile(f,np.int32).reshape((ngll,ngll,nel),order='F')

    # #read gll information
    filename = self.directory + '/gll_sem2d.tab'
    g = np.genfromtxt(filename)
    x, w, h = g[0,:], g[1,:], g[2:,:]
    self.mdict = {"nel" : nel,
            "npgeo" : npgeo,
            "ngnod" : ngnod,
            "npt" : npt,
            "ngll" : ngll,
            "coord" : coord,
            "ibool" : ibool,
            "x" : x,
            "w" : w,
            "h" : h,
            }
#

  def __read_header(self):
    """
    Read seismic header file of SEM2DPACK simulation.
    The method broadcasts the simulation parameters and
    receiver coordinates instances.
    """

    print ('Reading header file...')
    fname = self.directory + '/SeisHeader_sem2d.hdr'
    data = pd.read_csv(fname, header=0, nrows=1, dtype=str, delim_whitespace=True)
                      
                      
    self.dt = float(data['DT'].values[0])
    self.npts = int(data['NSAMP'].values[0])
    self.nsta = int(data['NSTA'].values[0])
    print (self.dt, self.npts, self.nsta)

    self.rcoord  = np.zeros( (self.nsta, 2) )
    with open(fname, 'r') as f:
      data = pd.read_csv(fname, names=('x','z'), delim_whitespace=True, header=2, nrows=self.nsta)
#       print (data)
      self.rcoord[:,0] = ( data['x'].values )
      self.rcoord[:,1] = ( data['z'].values )
     
      try: 
        line = pd.read_csv(fname, header=self.nsta+3, nrows=1, dtype=str, delim_whitespace=True)
        self.nsta_extra = int(line.columns[0])
        print ('Extra station number: ', self.nsta_extra)
        self.Elastic = False
        self.extra_coord  = np.zeros( (self.nsta_extra, 2) )
        data = pd.read_csv(fname, names=('x_extra','z_extra'), delim_whitespace=True, \
                           header=self.nsta+4, nrows=self.nsta_extra)
#         print (data)
        self.extra_coord[:,0] = ( data['x_extra'].values )
        self.extra_coord[:,1] = ( data['z_extra'].values )        
      except:
        print ('--- No extra station ---')
        pass
  ###


  @staticmethod
  def readField(fname):
    with open(fname,'rb') as f:
      field = np.fromfile(f,np.float32)
    return field
#

  def read_seismo(self,component='x'):
    if component == 'z': filename = self.directory + '/Uz_sem2d.dat'
    elif component == 'x': filename = self.directory + '/Ux_sem2d.dat'
    elif component == 'y': filename = self.directory + '/Uy_sem2d.dat'

    try :
      with open(filename, 'rb') as fid:
        veloc_array = np.fromfile(fid,np.float32)
    except : raise Exception('Velocity file does not exist')

    l = len(veloc_array)
    self.velocity = np.zeros((self.npts,self.nsta))

    if self.Dynamic: 
      limit=int(l/self.nsta)
      for i in range(limit):
        limit1 = i*self.nsta
        limit2 = (i+1)*self.nsta
        self.velocity[i,:] = veloc_array[limit1:limit2]
    else:
      limit = self.nsta
      for i in range(limit):
        limit1 = i*self.npts
        limit2 = (i+1)*self.npts
        self.velocity[:,i] = veloc_array[limit1:limit2]

    self.time = np.arange(self.velocity.shape[0])*self.dt
    return self.velocity
    ### 



  def read_seismos_rsf(self, ff=np.float32, jump=1):
    print ('Reading ', self.directory+'Seismos_rsf.dat')
    with open(self.directory+'Seismos_rsf.dat', 'rb') as fid:
      whole = np.fromfile(fid, ff) 
      data = whole.reshape(int(len(whole)/(self.nsta+1)), self.nsta+1)
      print ('data shape: ', data.shape)
      self.velocity = data[::jump, 1:]  
      # Elif: modif the code later to delete time array if unnecessary
      self.time = np.arange(self.velocity.shape[0])* self.dt* jump
    return
    ###



  def read_stress_strain(self, aktif=False):   
    ff = np.float32
    if self.Elastic or self.nsta_extra <= 0:
        print ('ERROR: Elastic conditions')  
        return
    if self.double_precision: ff = np.float64

    # Strain 
    filename = self.directory+ "/Strain_sem2d.dat"
    self.strain = read_binary(self, filename, ff)
    self.strain = 2.0* self.strain   # convert from epsilon to gamma

    # Strain 
    filename = self.directory+"/Stress_sem2d.dat"
    self.stress = read_binary(self, filename, ff)
    
    # Aktif: plasticity-surface number 
    if aktif:
        filename = self.directory+"/Nonlinear_surfaces_sem2d.dat"
        self.aktif = read_binary(self, filename, ff)   
        # print ('debug', self.aktif)

    # Max values
    self.max_strain = np.zeros(self.nsta_extra)
    if aktif: self.max_aktif = np.zeros(self.nsta_extra)
    for sta in np.arange(self.nsta_extra):
        self.max_strain [sta] = max( abs(self.strain[ :, sta]) )    
        if aktif: self.max_aktif [sta] = max( abs(self.aktif[ :, sta]) )    


  def read_nonlinear_backbone_data_tocorrect(self):

      filename = self.directory+ "/Nonlinear_soil_sem2d.dat"
      x = []; z = []; gref=[]; Gmod=[]
      with open (filename, 'r') as f:
          lines = f.readlines()
          nline = len(lines)
          for i in range(0, nline, 2):
              _x, _z, _gref = lines[i].split()[:]
              _Gmod = lines[i+1].split()[0]
              x.append(float(_x))
              z.append(float(_z))
              gref.append(float(_gref))
              Gmod.append(float(_Gmod))
      f.close()

      _x = np.array(x)
      _z = np.array(z)
      gref = np.array(gref)
      Gmod = np.array(Gmod)

      npts =_x.shape[0]
      self.nonlinear_curve = np.zeros((npts,4))
      self.nonlinear_curve[:,0] = _x[:]
      self.nonlinear_curve[:,1] = _z[:]
      self.nonlinear_curve[:,2] = gref[:]
      self.nonlinear_curve[:,3] = Gmod[:]
  ###
        
        
        
  ###
    


  @staticmethod
  def interp(field,coord):
    """
    Interpolates argument field over a meshgrid.
    Meshgrid size depends on argument coord.
    """
    xcoord = coord[:,0]
    zcoord = coord[:,1]
    nbx = len(xcoord)
    nbz = len(zcoord)
    ext = [min(xcoord), max(xcoord), min(zcoord), max(zcoord)]
    x,z=np.meshgrid(np.linspace(ext[0],ext[1],1000),np.linspace(ext[2],ext[3],1000),sparse=True)
    y = gd((xcoord,zcoord),field,(x,z),method='linear')
    y =np.flipud(y)
    return y
#

  @staticmethod
  def rinterp(x,y,z):
    # Set up a regular grid of interpolation points
    dy = y[1]-y[0] ;
    Nx, Ny = 2*max(x), max(y)/dy
    xi, yi = np.linspace(min(x), max(x), Nx), np.linspace(0, max(y), Ny)
    Xi,Yi  = np.meshgrid(xi,yi)
    zi = gd((x, y), z,( Xi, Yi), method='linear')
    return xi, yi, zi
#

  def plot_meshnode(self):
    filename = self.directory + '/MeshNodesCoord_sem2d.tab'
    nel = self.mdict["nel"]
    n
    with open(filename,'r') as f:
      nodes = np.genfromtxt(f)
#

  def plot_source(self,savefile=None,source_name=None,tap=0.025,fmax=10.0, tmax=2.0):
    from matplotlib import style
    style.use('ggplot')
    #if not isinstance(source_name,str):
    #  print('source file name must be str object')
    if source_name: source_name = source_name
    else: source_name = 'SourcesTime_sem2d.tab'

    source_file = self.directory + '/'+ source_name
    with open(source_file,'r') as src:
      amp = np.genfromtxt(src)
    # plot spectra
    dt = amp[1,0]-amp[0,0]
    spec,f = fourier(amp[:,1],dt,p=tap)

    fig = plt.figure(figsize=(8,6)); set_style(whitegrid=True, scale=1.0)
    # fig.subplots_adjust(wspace=0.3)
    ax1  = fig.add_subplot(121)
    ax2  = fig.add_subplot(122)
    ax1.plot(amp[:,0],amp[:,1])
    ax2.plot(f, spec)
    ax1.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax2.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax1.set_xlabel('Time [s]',fontsize=14) ; ax1.set_ylabel('velocity [$ms^{-1}$]',fontsize=14)
    ax2.set_xlabel('Frequency [Hz]',fontsize=14) ; ax2.set_ylabel('amplitude',fontsize=14)
    ax1.set_title('Source time function',fontsize=16)
    ax2.set_title('Source spectrum',fontsize=16)
    ax2.set_xlim(0.0,fmax)
    ax1.set_xlim(0.0,tmax)
    #plt.tight_layout
    if savefile : plt.savefig(savefile)
    plt.show()
#


  @staticmethod
  def plot_im(matrix,vmin,vmax,cmin,cmax,**kwargs):
    sns.set_style('whitegrid')
    fig = plt.figure(figsize=(8,6))
    ax.add_subplot(111)
    im = ax.imshow(matrix,cmap='jet',aspect='auto',interpolation='bilinear', \
                   vmin=vmin, vmax=vmax, origin='lower', extent=extent)
    if 'xlim' in kwargs   : ax.set_xlim(kwargs['xlim'][0], kwargs['xlim'][1])
    if 'ylim' in kwargs   :
      ax.set_ylim(kwargs['ylim'][0], kwargs['ylim'][1])
    else : ax.set_ylim(0.1,fmax)

    if 'ylabel' in kwargs : ax.set_ylabel(kwargs['ylabel'], fontsize=16)
    if 'xlabel' in kwargs : ax.set_xlabel(kwargs['xlabel'], fontsize=16)
    if 'title' in kwargs  : ax.set_title(kwargs['title'],fontsize=18)

    # colorbar
    cb = fig.colorbar(im, shrink=0.5, aspect=10, pad=0.01,\
                     ticks=np.linspace(cmin,cmax,cmax+1), \
                     boundaries=np.linspace(cmin,cmax,cmax+1))
    cb.set_label('Amplification', labelpad=20, y=0.5, rotation=90, fontsize=15)
    plt.show()
#


  def filter_seismo(self,fmin=0.0, fmax=10.0,ftype='lowpass', compo='x', isRSF=False):
    """
    filter seismograms.
    Inputs:
      -freqs[tuple][(0.1,10)] : corner frequencies of filter response
      -ftype[str][default=bandpass] : filter type
    Return:
      -Updates self.velocity array.
    """
    print ('Filtering seissmograms ...')
    if not self.velocity.size:
      if not isRSF: 
        self.read_seismo(component=compo)
      else:
        read_seismos_rsf()
    ##

    if ftype == 'lowpass':
      for sta in np.arange(self.nsta):
          self.velocity[:,sta] = lowpass(self.velocity[:,sta],fmax, df=1.0/self.dt)
    ##

    if ftype == 'highpass':
      for sta in np.arange(self.nsta):
          self.velocity[:,sta] = highpass(self.velocity[:,sta],fmin, df=1.0/self.dt)
    ##

  ###



  def plot_Vs(self,vs_br):
    from scipy.spatial.distance import pdist
    vsfile = self.directory + '/Cs_gll_sem2d.tab'
    with open(vsfile,'r') as v:
      vs_int = pd.read_csv(v,sep='\s+',names=['vs','x','z'])
    tmp = vs_int.drop_duplicates()
    self.vs_int = tmp.drop(tmp[tmp['vs']==vs_br].index)
    #dx = pdist(self.vs_int['x'][:,None],'cityblock')
    #dx = np.min(dx[np.nonzero(dx)])
    #dz = pdist(self.vs_int['z'][:,None],'cityblock')
    #dz = np.min(dz[np.nonzero(dz)])
    minx,maxx = np.min(self.vs_int['x']),np.max(self.vs_int['x'])
    minz,maxz = np.min(self.vs_int['z']),np.max(self.vs_int['z'])
    l = len(self.vs_int['x'])
    xi,zi = np.linspace(minx,maxx,l), np.linspace(minz,maxz,l)
    Xi,Zi = np.meshgrid(xi,zi)
    #plt.scatter(self.vs_int['x'],self.vs_int['z'],c=self.vs_int['vs'],cmap='jet')
    #plt.show()
    x = self.vs_int['x'].values
    z = self.vs_int['z'].values
    vs = self.vs_int['vs'].values
    y = gd((x,z),vs,(Xi,Zi),method='nearest')
    plt.figure()
    plt.imshow(y,cmap='jet',aspect='auto')
    plt.show()
    db.set_trace()
#

  def read_fault(self, ff=np.float32, LENTAG=1, is_rate_and_state=False):
      
      from distutils import util
      
      ''' Script to read FltXX files.
      Assuming that a single boundary output has been defined for the fault.
      to modify later for multiple fault boundaries...
      , also to modify for files with data > 5.'''

      BC = []; fault = {}
      for n, f in enumerate([self.directory+'/Flt'+('%02d' % i)+'_sem2d.hdr' for i in np.arange(1,7)]):
          found = os.path.exists(f)
          if found :
              BC.append(n+1)
              print ('Fault boundary: ', BC)
              break
      if not found : 
          print ('No Flt .hdr file found!')
          return
      
      
      elif not is_rate_and_state:     
          
          # Header file
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_sem2d.hdr'
          data = pd.read_csv(fname, names=('npts','ndat','nsamp','delta'), delim_whitespace=True, header=0, nrows=1)
          fault['npts'] = data['npts'].values[0]
          fault['ndat'] = data['ndat'].values[0]
          fault['nsamp'] = data['nsamp'].values[0]
          fault['delta'] = data['delta'].values[0]
          with open(fname, 'r') as f:
              line  = f.readlines()[2:3][0]
              fault['dat_names'] =  [el.replace('\n','').replace(' ','') for el in line.split(':')]
          data = pd.read_csv(fname, names=('x','z'), delim_whitespace=True, header=3)
          fault['x'] = data['x'].values
          fault['z'] = data['z'].values        

          # Init file
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_init_sem2d.tab'
          if os.path.exists(fname):
              data = pd.read_csv(fname, names=('st0','sn0','mu0'), delim_whitespace=True, header=None)
              fault['st0'] = data['st0'].values
              fault['sn0'] = data['sn0'].values
              fault['mu0'] = data['mu0'].values

          # Read fault data in a big matrix
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_sem2d.dat'
          if os.path.exists(fname):
              with open(fname, 'rb') as fid:
                  whole = np.fromfile(fid, ff) 
                  # BUG : nsamp is not correct inside the code !
                  try:
                      array = whole.reshape((2*LENTAG+fault['npts'], fault['ndat'], fault['nsamp']), order='F')
                  except:
                      fault['nsamp'] -=1
                      array = whole.reshape((2*LENTAG+fault['npts'], fault['ndat'], fault['nsamp']), order='F')

                  for j in np.arange(fault['ndat']):
                      print ('Assigning ', fault['dat_names'][j])
                      dat = fault['dat_names'][j]
                      data = array[LENTAG:LENTAG+fault['npts'], j, :]
                      fault [dat] = data
          fault['Time'] = np.linspace(0.0, fault['delta']* (fault['nsamp']), num=fault['nsamp'])

      
      elif is_rate_and_state:
          
          # Header file
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_sem2d.hdr'
          data = pd.read_csv(fname, names=('npts','ndat'), delim_whitespace=True, header=0, nrows=1)
          fault['npts'] = data['npts'].values[0]
          fault['ndat'] = data['ndat'].values[0]    
          with open(fname, 'r') as f:
              line  = f.readlines()[2:3][0]
              fault['dat_names'] =  [el.replace('\n','').replace(' ','') for el in line.split(':')]        
          data = pd.read_csv(fname, names=('x','z'), delim_whitespace=True, header=3)
          fault['x'] = data['x'].values
          fault['z'] = data['z'].values         
          
          
          # Init file
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_init_sem2d.tab'
          if os.path.exists(fname):
              fault['st0'] = np.genfromtxt(fname, usecols=0)
              fault['sn0'] = np.genfromtxt(fname, usecols=1)
              fault['mu0'] = np.genfromtxt(fname, usecols=2)
              fault['theta0'] = np.genfromtxt(fname, usecols=3)
              fault['V0'] = np.genfromtxt(fname, usecols=4)            
              fault['a'] = np.genfromtxt(fname, usecols=5)
              fault['b'] = np.genfromtxt(fname, usecols=6)
            
      
          # Data file I
          # Read fault data in a big matrix
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_sem2d.dat'
          if os.path.exists(fname):
              with open(fname, 'rb') as fid:
                  whole = np.fromfile(fid, ff) 
                  
                  nsamp = len(whole)/(2*LENTAG+fault['npts'])/fault['ndat']
                  print ('Guessed array size, nsamp: ', nsamp)
                  fault['nsamp'] = int(nsamp)
                  
                  array = whole.reshape((2*LENTAG+fault['npts'], fault['ndat'], fault['nsamp']), order='F')
                  
                  for j in np.arange(fault['ndat']):
                      print ('Assigning ', fault['dat_names'][j])
                      dat = fault['dat_names'][j]
                      data = array[LENTAG:LENTAG+fault['npts'], j, :]
                      fault[dat] = data                    
                  
                  
          # Data file II 
          # Read fault data in a big matrix
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_time_sem2d.tab'        
          if os.path.exists(fname):
              print(fname)
              fault['it'] = np.genfromtxt(fname, usecols=0)
              fault['dt'] = np.genfromtxt(fname, usecols=1)
              fault['t'] = np.genfromtxt(fname, usecols=2)
              fault['#EQ'] = np.genfromtxt(fname, usecols=3)
              dum = np.genfromtxt(fname, usecols=4, dtype=str)
              fault['isDyn'] = np.array( [bool(util.strtobool(d)) for d in dum] )             
              dum = np.genfromtxt(fname, usecols=5, dtype=str)
              fault['isSwitch'] = np.array( [bool(util.strtobool(d)) for d in dum] )
              dum = np.genfromtxt(fname, usecols=6, dtype=str)
              fault['isEq'] = np.array( [bool(util.strtobool(d)) for d in dum] )      

          # Data file III
          # Read potency and potency rate, currently only for out-of-plane (ndof=1)
          # for in-plane models change usecols and potency* arrays' size
          fname = self.directory+'/Flt'+str('%02d' % BC[0])+'_potency_sem2d.tab'        
          if os.path.exists(fname):
              print(fname)
              # POTENCY
              # out-of plane model -- compo 13
              # replace D by e for python 
              array = np.genfromtxt(fname,usecols=0, dtype=None,  encoding=None)
              array_fixed  = np.array( [float(dum.replace('D','e')) for dum in array])
              fault['potency'] = np.zeros((array_fixed.shape[0], 2))
              fault['potency_rate'] = np.zeros((array_fixed.shape[0], 2))    
              fault['potency'][:,0] = array_fixed
              
              # out-of plane model -- compo 23
              array = np.genfromtxt(fname,usecols=1, dtype=None,  encoding=None)
              array_fixed  = np.array( [float(dum.replace('D','e')) for dum in array])    
              fault['potency'][:,1] = array_fixed
              
              # POTENCY RATE
              # out-of plane model -- compo 13
              array = np.genfromtxt(fname,usecols=2, dtype=None,  encoding=None)
              array_fixed  = np.array( [float(dum.replace('D','e')) for dum in array])    
              fault['potency_rate'][:,0] = array_fixed
                  
              # out-of plane model -- compo 23
              array = np.genfromtxt(fname,usecols=3, dtype=None,  encoding=None)
              array_fixed  = np.array( [float(dum.replace('D','e')) for dum in array])    
              fault['potency_rate'][:,1] = array_fixed    
    
          
      self.fault = fault        

      # Find the earthquake
      cdt_beg = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], 1)==False)
      cdt_end = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], -1)==False)
      index = np.arange(0, len(self.fault['isDyn']))
      print ('Number of dynamic beginning and ending points:', len(index[cdt_beg]), len(index[cdt_end]))      
      return 
      
  ###



  def plot_cycles_time_step(self, savefig=False):

      plt.figure(figsize=(4,3))
      plt.xlabel('Simulation time (s)')    
      plt.ylabel('Time step dt (s)')      
      plt.semilogy(self.fault['t'], self.fault['dt'], 'k', marker='.')
      plt.grid()
      plt.tight_layout()
      if savefig: fig.savefig('/Users/elifo/Desktop/timestep.png', dpi=300)      
      plt.show()
  ###

  def plot_cycles_cumulative_slip(self,VW_halflen=9.5,savefig=False,jump_sta=1,jump_dyn=1,
                                      col_sta='gray',col_dyn='tomato'):

      size = min(self.fault['Slip'].shape[1], self.fault['isDyn'].shape[0] )
      slip = self.fault['Slip'][:, :size]
      isDyn = self.fault['isDyn'][:size]

      sta = slip[:, isDyn == False][:, ::jump_sta]
      dyn = slip[:, isDyn == True][:, ::jump_dyn]

      print ('WARNING: If weird output, try a smaller jump!')
      print ('shape of slip,sta,dyn: ', slip.shape,sta.shape,dyn.shape)

      fig = plt.figure(figsize=(8,4))
      for ii in enumerate(range(sta.shape[1])):
        plt.plot(sta[:, ii], self.fault['x'], c=col_sta,lw=0.1)
      #
      for ii in enumerate(range(dyn.shape[1])):
        plt.plot(dyn[:, ii], self.fault['x'], c=col_dyn,lw=0.1)
      #
      _xmin, _xmax = 0.0, max(slip[:,-1])
      plt.xlabel('Cumulated slip (m)')
      plt.ylabel('Along dip (m)')
      plt.hlines(VW_halflen, xmin=_xmin, xmax=_xmax, linestyle=':', color='g')
      plt.hlines(-VW_halflen, xmin=_xmin, xmax=_xmax, linestyle=':', color='g')
      if savefig: fig.savefig('/Users/elifo/Desktop/cumslip.png', dpi=300)  
      plt.show()    
  ###


  def plot_cycles_slip_rate(self, eq=1, is_normalisation=False, VW_halflen=9.5, _vmin=0.0, _vmax=2.0, 
                            savefig=False, VS_LVFZ=0.0, Lnuc=1.0, Vpl=6.34e-11, tmin=-1.0, tmax=-1.0, 
                            _cmap='rainbow'):

      import matplotlib.colors as colors
      
      # Find the earthquake
      cdt_beg = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], 1)==False)
      cdt_end = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], -1)==False)
      index = np.arange(0, len(self.fault['isDyn']))

      print ('Number of dynamic beginning and ending points:', len(index[cdt_beg]), len(index[cdt_end]))
      # choose the dynamic event range
      cdt1, cdt2 = index[cdt_beg][eq], index[cdt_end][eq]
      print ('cdt1, cdt2: ', cdt1, cdt2)

      # Slip rate and tim
      V = self.fault['Slip_Rate'][:, cdt1:cdt2] #npts, nt
      t = self.fault['t'][cdt1:cdt2]
      t -= t[0] # init time at zero

      # Shear stress: beginning vs end
      init = self.fault['Shear_Stress'][:,cdt1:cdt2][:,0] /1e6
      final = self.fault['Shear_Stress'][:,cdt1:cdt2][:,-1]/ 1e6

      init += self.fault['st0']/1e6
      final += self.fault['st0']/1e6
          
      # SLIP RATE
      data = V
      print ('Max slip rate: ', np.amax(V) )

      if is_normalisation: 
          data = V/ Vpl
          data = np.log10(data)
          # print ('After normal. max : ', np.amax(data) )
          xx = t* self.Vdamage/ Lnuc
          yy = self.fault['x']/Lnuc     

      else:
          xx = t
          yy = self.fault['x']/Lnuc
      #
                
          
      ### Plot
      print ('Plotting ...')
      fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]}, sharey=True)

      ###
      ax = axs[0]
      ax.set_title('Shear stress (MPa)')
      ax.set_ylabel('Along dip (Lnuc)')
      xmin, xmax = min(min(init), min(final)), max(max(init), max(final))
      ax.set_ylim(-2*VW_halflen/Lnuc, 2*VW_halflen/Lnuc)
      ax.hlines(y=VW_halflen/Lnuc, xmin=xmin, xmax=xmax, linestyle=':')
      ax.hlines(y=-VW_halflen/Lnuc, xmin=xmin, xmax=xmax, linestyle=':')
      ax.plot(init, yy, c='k', label='Initial stress')
      ax.plot(final, yy, c='red', label='Final stress')
      ax.grid()

      ###
      ax = axs[1]
      ax.set_title('Slip rate for event #'+ str(eq))
      xmin, xmax = tmin, tmax
      ext = [min(xx), max(xx), min(yy), max(yy)]      
      if tmax < 0.0: tmax = max(xx)
      if tmin < 0.0: tmin = min(xx)
      ax.set_xlim(tmin, tmax)
      ax.hlines(y=VW_halflen/Lnuc, xmin=xmin, xmax=xmax, linestyle=':')
      ax.hlines(y=-VW_halflen/Lnuc, xmin=xmin, xmax=xmax, linestyle=':')
      im = ax.imshow(data, extent=ext,
                 interpolation='nearest', cmap=_cmap, aspect='auto', vmin=_vmin, vmax=_vmax, origin='lower')
      ax.set_xlabel('t (V/Lnuc)')
      cb = fig.colorbar(im, orientation='vertical')
      cb.set_label('log(V/Vpl)')
      fig.set_size_inches(9.0, 4.5)
      plt.tight_layout()

      if savefig: fig.savefig('/Users/elifo/Desktop/event_'+str(int(eq))+'.png', dpi=300)        
      print('*')
  ###

  def get_static_iteration(self, eq=0):
      # Find the earthquake
      cdt_beg = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], 1)==False)
      cdt_end = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], -1)==False)
      index = np.arange(0, len(self.fault['isDyn']))
      print ('Number of dynamic beginning and ending points:', len(index[cdt_beg]), len(index[cdt_end]))
      # choose the dynamic event range
      if eq==0: 
          cdt1 = 0
          cdt2 = index[cdt_beg][eq]
      else:
          cdt1, cdt2 = index[cdt_end][eq-1], index[cdt_beg][eq]
      ###
      print ('Before eq', eq)
      print ('Between iterations: ', cdt1, cdt2)
      return cdt1, cdt2
  ###

  def get_dynamic_iteration(self, eq=0):
      # Find the earthquake
      cdt_beg = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], 1)==False)
      cdt_end = (self.fault['isDyn']==True) & (np.roll(self.fault['isDyn'], -1)==False)
      index = np.arange(0, len(self.fault['isDyn']))
      print ('Number of dynamic beginning and ending points:', len(index[cdt_beg]), len(index[cdt_end]))
      # choose the dynamic event range
      cdt1, cdt2 = index[cdt_beg][eq], index[cdt_end][eq]
      ###
      print ('Before eq', eq)
      print ('Between indices: ', cdt1, cdt2)
      return cdt1, cdt2
  ###
  
  def animate_cycles(self, cdt1, cdt2, jump=3, VW_halflen=1900.0, sleep=0.00001, it_cut=0,_alpha=0.1,_alphainc=0.1):
      import time
      import pylab as pl
      from IPython import display      
    
      Vf = self.fault['Slip_Rate'][:, cdt1:cdt2]    
      Vf_max = np.array( [max(abs(v))  for v in Vf.T] )      
      slip = self.fault['Slip'][:, cdt1:cdt2]      
      it = self.fault['it'][cdt1:cdt2]
      dt = self.fault['dt'][cdt1:cdt2]
      x = self.fault['x']
      dum = self.fault['Shear_Stress'][:, cdt1:cdt2]
      stress = dum+ self.fault['st0'][:, None]
      isDyn = self.fault['isDyn'][cdt1:cdt2]    

      if it_cut==0: it_cut=max(it)
      print ('Max it: ', it_cut)  
      
      fig = plt.figure(figsize=(12, 8))    

      ax1 = plt.subplot(131)
      plt.title('Stress (MPa)')
      plt.grid()
      ax1.set_xlim(50.0, 120.0)
      ax1.set_ylim(-2*VW_halflen, 2*VW_halflen)
      ax1.axhline(y=VW_halflen, linestyle=':')
      ax1.axhline(y=-VW_halflen, linestyle=':')                
      
      # subplot: slip rate
      ax2 = plt.subplot(132)
      plt.grid()
      plt.title('Slip rate (m/s)')
      ax2.set_xlim(1e-20, 15.0)
      ax2.set_ylim(-2*VW_halflen, 2*VW_halflen)
      ax2.axhline(y=VW_halflen, linestyle=':')
      ax2.axhline(y=-VW_halflen, linestyle=':')
            
      # subplot: slip rate
      ax3 = plt.subplot(133)
      plt.grid()
      plt.title('Slip (m)')
      ax3.set_ylim(-2*VW_halflen, 2*VW_halflen)
      ax3.axhline(y=VW_halflen, linestyle=':')
      ax3.axhline(y=-VW_halflen, linestyle=':')         
            
      ii=0
      for v,d, _it,_dt, _stress, _isDyn in zip(Vf.T[::jump], slip.T[::jump], it[::jump], dt[::jump], stress.T[::jump], isDyn):
          ii += 1
          print ('Step and index: ', _it, ii)
          ax2.semilogx(v, x, 'peru', alpha=min(_alpha+_alphainc*ii, 1.0))
          ax1.plot(_stress/1e6, x,   alpha=min(_alpha+_alphainc*ii, 1.0), c='k')
          ax3.plot(d, x, 'brown', alpha=min(_alpha+_alphainc*ii, 1.0))
          display.clear_output(wait=True)          
          # re-drawing the figure
          fig.canvas.draw()
          # to flush the GUI events
          fig.canvas.flush_events()        
          display.display(pl.gcf())
          time.sleep(sleep)          
          if _it >=it_cut: break
      ###
  #


  def plot_2D_slip_rate(self,  save=False, figname='2d_fault', cmap='magma', **kwargs):
    ''' Spatio-temporal plot for slip rate along the fault line.
    Only positive x stations are used.
    '''

    print ('Plotting the spatiotemporal graph of slip rate...')
    fig = plt.figure(figsize=(6,6)); set_style(whitegrid=False, scale=1.0)
    
    ax = fig.add_subplot(111)
    ax.set_xlabel('Time ($L_{c}/V_{s}^{host}$)')
    ax.set_ylabel('Distance ($L_{c}$)')
    if 'ylimits' in kwargs: 
      ylimits = kwargs['ylimits']
      ax.set_ylim(ylimits[0], ylimits[1])


    # Data mesh - x
    time = self.fault['Time']; x = time
    # Data mesh - y
    xcoord = self.fault['x']
    jj = np.ravel(np.where(xcoord >= 0.0))
    y = xcoord[jj] ;  ylim = [0.0, max(y)]
    # if 'ylim' in kwargs: ylim = kwargs['ylim']
    ext = [min(x), max(x), ylim[0], ylim[1]]   
    # Data
    z = self.fault['Slip_Rate'][jj,:]  
    index = np.where( (y >= ylim[0]) & (y <= ylim[1]) )[0]
    vmin = 1.e-2; vmax = z[index].max() # for LogNormal scale
    if 'vmax' in kwargs: vmax = kwargs['vmax']    
    print ('Min and Max of data: ', z.min(), z.max())
    print ('Min and Max of chosen domain: ', z[index].min(), z[index].max())
    z [z < vmin] = vmin

    # tit = 'Max = '+ str('%.2f' %  z[index].max())
    # tit += ' at distance = '+ str('%.2f' %  y[np.where(z == z[index].max()) [0] ] )
    # ax.set_title(tit)
    im = ax.imshow(z, extent=ext, cmap= cmap, origin='lower',\
                 norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='bicubic')

    c = plt.colorbar(im, fraction=0.046, pad=0.1, shrink=0.3)
    c.set_clim(vmin, vmax); c.set_label('Slip rate ($V_{dyn}$)')

    plt.tight_layout()

    if save: plt.savefig(figname+'.pdf', dpi=300)
    plt.show(); plt.close()
#


  def plot_slip_rate(self, dist=1., save=False, figname='fault_data'):

    print ('Plotting the outputs at a distance of ', dist )
    fig = plt.figure(figsize=(8,6))
    sns.set_style('whitegrid')
    plt.subplots_adjust(top=0.88, bottom=0.11, left=0.09, right=0.96,\
                  hspace=0.3, wspace=0.3)

    time = self.fault['Time']
    xcoord = self.fault['x']
    jj = np.ravel( np.where(abs(xcoord-dist) < 1.e-5 ) ) [0]
    # print ('Found index : ', jj, xcoord[jj])

    slip = self.fault['Slip'][jj,:]
    srate = self.fault['Slip_Rate'][jj,:]
    mu = self.fault['Friction'][jj,:]
    tau = self.fault['Shear_Stress'][jj,:]
    sigma = self.fault['Normal_Stress'][jj,:]
    sigma_0 = abs(self.fault['sn0'][jj])
    tau_0 = self.fault['st0'][jj]


    # Slip-time function
    ax = fig.add_subplot(221)
    plt.suptitle('Fault data at distance ($L_{c}$): '+ str('%.3f' % dist), fontsize=18 )
    ax.set_ylabel('Slip ()',fontsize=16)
    ax.plot(time, slip, color='black')

    # Slip-rate function
    ax = fig.add_subplot(223)
    ax.set_xlabel('Time ()',fontsize=16)
    ax.set_ylabel('Slip rate ()',fontsize=16)
    ax.plot(time, srate, color='black')

    # Friction coefficient mu
    ax = fig.add_subplot(222)
    ax.set_ylabel('Friction coefficient',fontsize=16)
    ax.plot(time, mu, color='black')

    # Stress and strength changes
    ax = fig.add_subplot(224)
    ax.set_xlabel('Time ()',fontsize=16)
    # ax.plot(time, sigma, color='black', label='Normal stress')
    # ax.plot(time, tau, color='red', label='Shear stress')
    ax.plot(time, sigma+ sigma_0, color='black', label='Normal stress')
    ax.plot(time, tau+tau_0, color='red', label='Shear stress')    
    ax.plot(time, mu*sigma_0, color='blue', label='Shear strength')
    ax.legend(loc='best')
    if save: plt.savefig(figname+'.png', dpi=300)
    plt.show(); plt.close()
#



  def plot_snapshot_tests(self,fname,interval, vmin=-1.e-10, vmax=1.e-10, save=False,outdir='./',
            show=False, nsample=1000, cmap='seismic',\
            ylabel='Width / $L_{c}$', xlabel='Length / $L_{c}$', \
            lims=None, switch_coord=False, logscale=False, normaliseby=1.0,roll_xcoord=False):

    print ('Plotting snapshots...')
    filename = self.directory+ fname
    print ('reading ...', filename)
    field = self.readField(filename)

    coord = self.mdict["coord"]
    xcoord, zcoord = coord[:,0]/normaliseby, coord[:,1]/normaliseby
    # switch coord for anti-plane cycle simulations 
    # using 'x' along fault dip
    if switch_coord: 
      xcoord, zcoord = coord[:,1]/normaliseby, coord[:,0]/normaliseby
    if roll_xcoord:
      xcoord = max(xcoord)- xcoord
    nbx = len(xcoord)/4 ; nbz = len(zcoord)/4
    ext = [min(xcoord), max(xcoord), min(zcoord), max(zcoord)]
    print ('Model extent: ', ext)

    # Test for Nepal simulations
    # nearest is much faster than linear !
    x,z = np.meshgrid(np.linspace(ext[0],ext[1],nsample),np.linspace(ext[2],ext[3],nsample))
    y = gd((xcoord,zcoord),field,(x,z),method='nearest')
    print ('Min, Max of Field: ', np.amin(field), np.amax(field))
    print('flipud ...')
    y = np.flipud(y)

    #
    fig = plt.figure()
    sns.set_style('whitegrid')
    ax = fig.add_subplot(111)
    if lims != None:
      ax.set(xlim=(lims[0],lims[1]), ylim=(lims[2],lims[3]))

    # Fault rupture outputs
    if not logscale:
      im = ax.imshow(y, extent=[min(xcoord), max(xcoord), min(zcoord), max(zcoord)], \
        vmin=vmin, vmax=vmax, cmap=cmap, aspect='auto')
    else:
      im = ax.imshow(np.log10(abs(y)), extent=[min(xcoord), max(xcoord), min(zcoord), max(zcoord)], \
        vmin=vmin, vmax=vmax, cmap=cmap, aspect='auto')
      print ('Min, Max of log10(field): ', np.amin(np.log10(abs(y))), np.amax(np.log10(abs(y))) )

    # custom cycle plots
    # add VW limits
    plt.axhline(y=self.VW_halflen/normaliseby, c='snow',linestyle=':')
    plt.axhline(y=-self.VW_halflen/normaliseby, c='snow',linestyle=':')
    try:
      plt.axvline(x=self.LVFZ/normaliseby, c='k',linestyle=':')
    except:
      pass

    # Adding rectangle to restrain the fault area (optional)
    # ax.add_patch(Rectangle((-15.0, -1.5),30., 3.,alpha=1,linewidth=1,edgecolor='k',facecolor='none'))

    plt.ylabel(ylabel); plt.xlabel(xlabel)
    c = plt.colorbar(im, fraction=0.046, pad=0.1,shrink=0.4)
    # c.set_clim(vmin, vmax)
    c.set_label('Amplitude')
    tit = 'Snapshot at t (s)= '+ str(interval)
    if interval < 0: tit = 'File: '+ fname    
    if logscale:
      tit += '   Min- Max vel. amplitude = '+ str('%.2f' % (np.amin(np.log10(abs(y)))))
      tit += '   '+ str('%.2f' % (np.amax(np.log10(abs(y)))))
    else:
      tit += '   Min- Max vel. amplitude = '+ str('%.2f' % (min(abs(field))))
      tit += '   '+ str('%.2f' % (max(abs(field))))      

    plt.title(tit); plt.tight_layout()
    if save: plt.savefig(fname+'.png',dpi=300) # save into current directory
    if show: plt.show()
    plt.close()
    print ('*')
#


  def animate_fault(self, compo='x', field='v', t_total=10.01, itd=500,\
                      vmin=-2.5, vmax=2.5, ready=False, digit=2,cmap='seismic',\
                      ibeg=0, iend=1, interval=-1,jump=1,xlabel='',ylabel='',\
                      lims=None,switch_coord=False,logscale=False, \
                      normaliseby=1.0,roll_xcoord=False):
    ''' Preparing snapshots and their gif in the current path. '''
    # Make snapshots from binary files

    if iend < 0:
      interval = round(self.dt, digit)* float(itd)
      iend = int(t_total/interval)+ 1
      print ('*** interval and total number: ', interval, total)
    if interval <0 :
      print ('Provide interval: interval*i is time of snapshots!')
      print ('*')

    files = []
    for i in range(ibeg, iend, jump):
      n = str('%03d' % i)
      fname = field+compo+'_'+n+'_sem2d.dat'
      if not ready:
        # self.plot_snapshot(fname, interval*i, vmin=vmin, vmax=vmax, save=True, show=False) 
        self.plot_snapshot_tests(fname, interval*i, vmin=vmin, vmax=vmax, save=True, show=False, cmap=cmap,\
                                  xlabel=xlabel, ylabel=ylabel,lims=lims,\
                                  switch_coord=switch_coord,logscale=logscale,\
                                  normaliseby=normaliseby,roll_xcoord=roll_xcoord)     

      files.append(fname+'.png')
    # Animate the plot_snapshot_testspshots
    create_gif(files, 1.5)
#   



  def plot_fronts(self, Dc=1.0, Veps=1.e-10, head=True, tail=True, diff=False,\
                  eps=1.e-3, d_elem=0.25, sigma_gaussian=3, \
                  save=False, fname='fronts' , debug=False, **kwargs):
 
    ''' Only for the grid points btw x=0 and x=xmax '''

    slip = self.fault['Slip']
    srate = self.fault['Slip_Rate']
    time = self.fault['Time']
    xcoord = self.fault['x']

    # Maximum distance to be used 
    xmax = xcoord.max() 
    if 'xmax' in kwargs: xmax = kwargs['xmax']

    # Compute passages of rupture front and tail
    dt = time[1]-time[0]
    nn = np.where( (xcoord >=0.0) & (xcoord <= xmax))[0]
    xx = xcoord[nn]; sr = srate[nn,:]; sl = slip[nn,:]
    self.fault['x_rupt'] = xx

    # rupture front (here long way to avoid noise !!!)
    ii=[]; Trupt = np.zeros( len(xx)); t_limit= np.zeros( len(xx))
    for n in np.arange(len(xx)):
      Vbeta = 1.e-10
      t_peak, multipulse = find_tpeak(time, sr[n,:], Vbeta)
      if multipulse and debug: 
         print ('*** MULTIPLE PULSE at (x,t): ', xx[n], t_peak)
         print ('Define a smaller xmax (by default, it equals model length')
      iii = np.where ( (time < t_peak.min()) & (sr[n,:] < Vbeta) )[0]
      # iii = np.where ( (time < t_peak) & (sr[n,:] < Vbeta) )[0]
      if len(iii) > 0 :
        index = iii[::-1][0]
        t = time [index]    
        t_limit[n] = t
        Trupt[n] = t+ dt* (Veps-sr[n,index])/ (sr[n,index+1] - sr[n,index])
        #
        if not multipulse and  t <= t_limit[n-1]  :
          # this part is to avoid the numerical prb
          # where zero value > Vbeta !
          while n > 0 and Vbeta < 0.15 : #1.e-2 :
             if debug : print ('PRECISION problem at (x,t) : ', xx[n], t, Vbeta )
             Vbeta *= 2.0 
             iii = np.where ( (time < t_peak) & (sr[n,:] < Vbeta) )[0]
             index = iii[-1]
             t = time [index] 
             Trupt[n] = t+ dt* (Veps- sr[n,index])/ (sr[n,index+1] - sr[n,index])
        ###
      else:
        print ('Rupture front does not reach to distance ', xx[n])
        break 
    #

    # end of process zone (tail)     
    kk = np.array ( [ np.where( sl[n,:] < Dc )   for n in np.arange(len(xx)) ] )
    n=0; idxx=0; Tproz = np.zeros( len(kk) )
    while n < len(kk):
      for k in kk:
        if len (k[0])>0 :
          idxx = k[0][-1] 
          if k[0][-1] < len(time)-1:
            Tproz[n] = time[idxx]
            Tproz[n] += dt* (Dc - sl[n,idxx])/ (sl[n,idxx]-sl[n,idxx-1])
        n+=1
    #
    self.fault['Tproz'] = Tproz
    self.fault['Trupt'] = Trupt  
    # Difference
    jj = np.where (abs(xx) < eps)[0]
    diff0 = Tproz[jj] - Trupt[jj]


    # get rupture speed by derivative
    kk = np.where( (Trupt  > 0.0) ) [0]
    V = np.diff(xx[kk])/ np.diff(Trupt[kk])
    # smooth rupture speed
    Vrupt = scipy.ndimage.filters.gaussian_filter1d(V, sigma_gaussian)
    self.fault['Vrupt'] = Vrupt
    # self.fault['Vdum'] = V
    

    # Critical distance computation
    fnc = scipy.interpolate.griddata(Tproz, xx, \
              Trupt[np.where( (Trupt > Tproz[0]) & (Trupt < Tproz.max()) )], method='linear')
    xxx = xx [ np.where( (Trupt > Tproz[0]) & (Trupt < Tproz.max())  ) ]    
    Lc  = xxx - fnc
    Lplot = True
    self.fault['Lc'] = Lc




    # Plot
    print ('Plotting the rupture fronts...')
    print ('Initial difference btw rupture front and tail (s): ', diff0 )
    fig = plt.figure (figsize=(12,6))
    sns.set_style('whitegrid')
    plt.subplots_adjust(left=0.085, right=0.955, wspace=0.285, hspace=0.205)
    #
    ax = fig.add_subplot(131)
    ax.set_title('Initial pulse width (in time): '+ str('%.2f' % diff0), fontsize=16)
    ax.set_xlabel('Distance along the fault ()', fontsize=16)
    ax.set_ylabel('Time ()', fontsize=16)
    ax.set_xlim([0.0, xmax]); ax.set_ylim([0.0, 1.1* Tproz.max()])
    if 'xlim' in kwargs:  ax.set_xlim([kwargs['xlim'][0], kwargs['xlim'][1]])
    if 'tlim' in kwargs:  ax.set_ylim([kwargs['tlim'][0], kwargs['tlim'][1]])
    if head : ax.plot(xx, Trupt, color='black', label= 'Rupture front');
    if tail : ax.plot(xx, Tproz, color='royalblue', label='Tail of process zone'); 
    if diff : ax.plot(xx, Tproz-Trupt, color='black', label='Difference', linestyle='--', alpha=0.7);   
    plt.legend(loc='best')
    #
    if Lplot:
      ax = fig.add_subplot(132)
      ax.set_title('Minimum process zone: '+ str('%.2f' % min(Lc)), fontsize=16)
      ax.set_xlim(0.0, xmax);  ax.set_ylim(0.0, Lc.max())
      if 'xlim' in kwargs:  ax.set_xlim([kwargs['xlim'][0], kwargs['xlim'][1]])
      ax.set_xlabel('Distance along the fault ()', fontsize=16)
      ax.set_ylabel('Process zone', fontsize=16)
      ax.plot(xxx, Lc, color='k')
      labo = 'h = '+ str(d_elem)
      plt.axhline(y=d_elem, xmin=0.0, xmax=xmax, linestyle='--', color='gray',lw=5., label=labo)
    #
    ax = fig.add_subplot(133)
    ax.set_xlim([0.0, xmax]); ax.set_ylim(0.0, 1.1*Vrupt.max())
    ax.set_title('Max. speed: '+ str('%.2f' % Vrupt.max()), fontsize=16)   
    if 'xlim' in kwargs:  ax.set_xlim([kwargs['xlim'][0], kwargs['xlim'][1]])  
    ax.set_xlabel('Distance along the fault ()', fontsize=16)
    ax.set_ylabel('Rupture speed ()', fontsize=16)    
    
    ax.plot(xx[kk][:-1], abs(Vrupt), color='k')
    ax.plot(xx[kk][:-1], abs(V), color='gray', alpha=0.5)

    # ax.plot(xx[kk][:-1], Vrupt, color='k')
    # ax.plot(xx[kk][:-1], V, color='gray', alpha=0.5)
    if save: plt.savefig(fname+'.png',dpi=300) # save into current directory
    plt.show(); plt.close()
###




###