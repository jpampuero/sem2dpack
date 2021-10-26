# Script that has subroutines to compute FFT, smooth, etc

import numpy as np
# import ipdb as db
import sys

# Tapering with a Hanning window

def taper(x,p):
   if p <= 0.0:
      return x
   else:
      f0 = 0.5
      f1 = 0.5
      n  = len(x)
      if (p*n) < 1 : return x
      nw = int(p*n)
      ow = np.pi/float(p*n)

      w = np.ones( n )
      for i in range( nw ):
         w[i] = f0 - f1 * np.cos(ow*i)

      for i in range( n-nw,n ):
         w[i] = 1.0 - w[i-n+nw]

      return x * w

# Computing the fft of the data

def fourier(x,dt,p=0.0,remove_mean=False):

   # Removing the mean
   if remove_mean:
     x = x - np.mean(x)

   # Tapering the data
   if p > 0.0: x = taper(x,p)

   # FFT
   s  = np.abs( dt * np.fft.fft(x) )
   df = 1.0 / ((len(x)-1) * dt)
   nf = int( np.ceil( len(x)/2.0 ) ) + 1
   f  = np.arange( nf ) * df

   return s[:nf],f[:nf]

# Konno-Ohmachi smoothing
# Smoothing is done between 0 and fmax
# If fmax > Nyquist, all FFT is smoothed
# Too long in Python!!!

def ko(y,dt,dx,fmax=10.0,bexp=40.0):

  # dx is delta_frequency here
  if fmax > 1.0/(2.0 * dt):
      nx = len(y)
  else:
      nx = int(fmax/dx)

  fratio = 10.0**(2.5/bexp)
  ys     = np.zeros( nx )
  ys[0]  = y[0]

  for ix in range( 1,nx ):
     fc  = float(ix)*dx
     fc1 = fc/fratio
     fc2 = fc*fratio
     ix1 = int(fc1/dx)
     ix2 = int(fc2/dx) + 1
     if ix1 <= 0:  ix1 = 1
     if ix2 >= nx: ix2 = nx
     a1 = 0.0
     a2 = 0.0
     for j in range( ix1,ix2 ):
        if j != ix:
           c1 = bexp * np.log10(float(j)*dx/fc)
           c1 = (np.sin(c1)/c1)**4
           a2 = a2+c1
           a1 = a1+c1*y[j]
        else:
           a2 = a2+1.0
           a1 = a1+y[ix]
     ys[ix] = a1 / a2

     fs = np.arange(nx) * dx

  return fs, ys

def ko2(raw_signal,freq_array,smooth_coeff=40,progress_bar=False):
  x = raw_signal  # shorten variable names...
  f = freq_array
  f_shifted = f/(1+1e-4)  # shift slightly to avoid numerical errors
  b = float(smooth_coeff)

  if len(x) != len(f):
    print('Length of input signal and frequency array must be the same.')
    sys.exit()

  L = len(x)
  y = np.zeros(L)  # pre-allocation of smoothed signal

  if progress_bar == True:
    progress_bar_width = 40  # width of progress bar (40 characters)
    print('\n|------------    Progress    ------------|') # reference bar
    sys.stdout.write('|')

  # =======  Moving window smoothing: fc from f[1] to f[-2]  ========
  for i in range(L):
    if progress_bar and (np.remainder(i,L/progress_bar_width) == 0):
      sys.stdout.write('|')  # prints "|" without spaces or new lines

    if (i == 0) or (i == L-1):
      continue  # skip first and last indices for now

    fc = f[i]  # central frequency
    w = np.zeros(L)  # pre-allocation of smoothing window "w"

    z = f_shifted / fc  # "z" = dimensionless frequency, normalized by fc
    w = (np.sin(b * np.log10(z)) / b / np.log10(z)) ** 4.0
    w[np.isnan(w)] = 0  # replace NaN's with 0

    y[i] = np.dot(w,x) / np.sum(w)  # apply smoothing filter to "x"

  y[0] = y[1]  # calculate first and last indices
  y[-1] = y[-2]

  if progress_bar:
    sys.stdout.write('|\n')

  return y

def rtrend(x):
    ind = np.arange( len(x) )
    r   = np.polyfit(ind,x,1)
    fit = np.polyval(r,ind)
    x   = x - fit

    return x

# FFT Ratio

def fft_ratio(x,y,dt,smooth=40,taper=0.05,fmax=10):
    f, ax = fourier(x,dt,taper)
    f, ay = fourier(y,dt,taper)

    if smooth > 0:
        fs, ax = ko(ax,dt,f[1]-f[0],smooth,fmax)
        fs, ay = ko(ay,dt,f[1]-f[0],smooth,fmax)
        return fs, ay / ax
    else:
        return f, ay / ax

# Bitwise version

def next_power_of_2(n):
    """
    Return next power of 2 greater than or equal to n
    """
    return 2**(n-1).bit_length()

