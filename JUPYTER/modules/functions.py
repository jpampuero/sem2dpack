"""
Module containing functions necessary for post processing results
obtained from SEM2D simulation.

@author : Flomin TCHAWE

"""

import numpy as np
from matplotlib.mlab import griddata
from houches_fb import *
from filters import *
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import ipdb as db
import os
import multiprocessing as mp
#from __future__ import division


def read_header(filename,extra=False):
   f = open(filename, 'r')
   f.readline()
   string = f.readline()
   word   = string.rstrip(" ").split()	

   dt    = float(word[0])
   niter = int(word[1])
   nstat = int(word[2])
    
   # Seismos 
   f.readline()
   coord  = np.zeros((nstat,2))
   for stat in np.arange(nstat):
      string = f.readline()
      word   = string.rstrip(" ").split()			
      # x-coord
      coord[stat,0] = float(word[0])
      # z-coord
      coord[stat,1] = float(word[1])
    
   #extra station
   if extra :
      nxsta = f.readline()
      nxsta = int(nxsta)
      f.readline()
      xstat_coord = np.zeros((nxsta,2))
    
      for sta in range(nxsta):
         xtra = f.readline()
         word = xtra.rstrip(" ").split()
         xstat_coord[sta,0] = float(word[0])
         xstat_coord[sta,1] = float(word[0])
   else:
       xstat_coord=None
   
   return dt, niter, nstat, coord, xstat_coord
    

def read_seismo_dynamic(filename,direction='x'):
   if filename[-1] != '/': filename += '/'
   if direction == 'z': veloc_file = filename + "Uz_sem2d.dat"
   else: 
      print('default x direction is used')
      veloc_file = filename + "Ux_sem2d.dat"
          
   header_file = filename + "SeisHeader_sem2d.hdr"

   try:
      with open(veloc_file, 'rb') as fid:
         veloc_array = np.fromfile(fid, np.float32)   
         # veloc_array = np.fromfile(fid, np.float64) 
   except: raise Exception('Velocity file does not exist')

   dt,niter,nsta,coord = read_header(header_file)[0:4]

   l = len(veloc_array)
   nstat = int(l/niter)
   veloc_stats = np.zeros((niter,nstat))
   acc = np.zeros((niter,nstat))
   for i in range(int(len(veloc_array)/nstat)):
      limit1 = i*nstat
      limit2 = (i+1)*nstat
      veloc_stats[i,:] = veloc_array[limit1:limit2]
   #---
   # time vector
   t = np.arange(veloc_stats.shape[0])*dt

   return veloc_stats,t,coord
    

def tf(signal, dt, blim, xcoord, nsurf_recv,fmax,proc=1,smooth=False,*argv):
  """
    computes the transfer function of a sedimentary basin.
    
    ** param **
      - signal   : (2d array) recievers seismograms, contains bedrock
	         and basin recievers
      - dt       : (float) signals time step
      - blim : (tuple) x-coordinates of basin-bedrock limits
      - xcoord   : (1d-array) x-coordinates of recievers signal
      - nsurf_recv : (int) number of surface recievers
      - fmax     : (float) maximum frequency for representation
      - proc     : (int)  number of processors for parallel
	           implementation
      - smooth   : (bool) smoothing option, if True konno-ohmachi
	           smoothing function is applied to signals

  """

  #----------- 
  # Check input paramters
  
  if signal.ndim != 2 :
    print("Input array must be a 2-dimensional array")
  
  if (blim[0] < min(xcoord)) or (blim[1] > max(xcoord)) :
    print("basin-bedrock limits not correctly provided")
  
  x_ind = np.where(xcoord<blim[0])[0]
  y_ind = np.where(xcoord>blim[1])[0]
  rock_ids = np.append(x_ind,y_ind)
  nt, nx = signal.shape

  #------------
  #set run type option (sequential or parallel)

  if proc < 2 :    # sequential
    # compute bed rock fft
    rock_signal = []
    for i in rock_ids:
      rock,f = fourier(signal[:,i],dt,0.025)
      rock_signal.append(rock)
    rock_signal = np.array(rock_signal)
    rockm = np.mean(rock_signal,axis=0) 

    # compute 2D spectral ratio
    basin = []
    for i in range(nx):
      amp,f = fourier( signal[:,i], dt, 0.025)
      basin.append( amp / rockm )
    
  else :
    pool = mp.Pool(processes=proc)
    results = [pool.apply_async(fourier,args=(signal[:,i],dt,0.025)) for i \
               in rock_ids]
    mrock = [p.get()[0] for p in results]
    f = results[0].get()[1]
    del results
    
    mrock = np.array(mrock)
    rockm = np.mean(mrock,axis=0)

    # Smoothing the spectrum
    df    = f[1] - f[0]
    if smooth : fs, rockm = ko(rockm,dt,df,fmax)
    
    # Computing 2D spectral ratio 

    results = [pool.apply_async(fourier,args=[signal[:,i],dt,0.025]) for i \
              in range(nx)]
    basin = [p.get()[0]/rockm for p in results]
    
  basin = np.array(basin)

  if argv:
    ftfig=plt.figure()
    ax=ftfig.add_subplot(111)
    #ax.plot(f,basin[2],'r')
    ax.plot(f,rockm)
    plt.ion()
    plt.show() 
    db.set_trace()
    
  # interpolating TF 2D to a uniform grid
  x = []
  y = []
  z = []
  for i in range(nsurf_recv):
    for j in range( len(f) ):
      x.append( xcoord[i] )
      y.append( f[j] )
      z.append( basin[i,j] )
    
  z = np.array(z)
  xi, yi, zi = interp(x, y, z)
  return zi,basin, f

def interp(x,y,z):
  # Set up a regular grid of interpolation points
  dy = y[1]-y[0]
  Nx, Ny = 2*max(x), max(y)/dy 
  xi, yi = np.linspace(0, max(x), Nx), np.linspace(0, max(y), Ny)

  zi = griddata(x, y, z, xi, yi, interp='linear')

  return xi, yi , zi

def make_colormap(seq):
  """Return a LinearSegmentedColormap
  seq: a sequence of floats and RGB-tuples. The floats should be increasing
  and in the interval (0,1).
  """
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
    
    
 
def readStrainStressparam(filename,niter,nstat):
   with open(filename,'rb') as f:
      array = np.fromfile(f,np.float32)
   l = len(array)
   data = np.zeros((nstat,niter))
    
   for i in range(l/nstat):
      infli = i*nstat
      supli = (i+1)*nstat
      data[:,i]=array[infli:supli]
	       
   return data
   

def divergence(f):
   """
   The divergence of a vector field is :
      divergence(f) = dfx/dx + dfy/dy + dfz/dz + ...
   """
   num_dims = len(f)
   return np.ufunc.reduce(np.add, [np.gradient(f[i], axis=i) for i in range(num_dims)])


def correlate(trace1, trace2, t, maxlag=100, plot=False):
  """
  This function computes the correlation function of trace1 and trace2 as a function of time.
  The maximum time shift, maxlag, is the maximum number of index values by which the two discrete time series
  are shifted with respect to each other. If we only want to determine differential traveltimes, maxlag
  can be chosen pretty small.
  """
   
  #- Initialisations. ------------------------------------------
  time_index=np.arange(-maxlag,maxlag+1)
  tcc=time_index*(t[1]-t[0])
  cc=np.zeros(len(tcc))

  nt=len(t)
     
  #- Compute correlation function. -----------------------------
  for k in time_index:
    if k>0:
      cc[k+maxlag]=np.sum(trace1[k:nt]*trace2[0:nt-k])
    else:
      cc[k+maxlag]=np.sum(trace1[0:nt+k]*trace2[-k:nt])
                  
  #- Plot if wanted. -------------------------------------------
  if plot==True:
    plt.plot(tcc,cc)
    plt.show()

  return cc, tcc

