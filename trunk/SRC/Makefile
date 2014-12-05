# SEM2DPACK version 2.3.8 -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
#                            with emphasis on computational seismology and earthquake source dynamics.
# 
# Copyright (C) 2003-2007 Jean-Paul Ampuero
# All Rights Reserved
# 
# Jean-Paul Ampuero
# 
# California Institute of Technology
# Seismological Laboratory
# 1200 E. California Blvd., MC 252-21 
# Pasadena, CA 91125-2100, USA
# 
# ampuero@gps.caltech.edu
# Phone: (626) 395-6958
# Fax  : (626) 564-0715
# 
# http://web.gps.caltech.edu/~ampuero/
# 
# This software is freely available for academic research purposes. 
# If you use this software in writing scientific papers include proper 
# attributions to its author, Jean-Paul Ampuero.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# 
#== BEGIN USER SETTINGS =============================================

# 1. Modify next line if you need to place the executable elsewhere:
EXEC = $(HOME)/bin/sem2dsolve

# 2. Uncomment/add/modify the lines corresponding to your compiler:

# GNU 
#F90 = g95
#OPT = 
#F90 = gfortran
#OPT = -O3 -Wall
#OPT = -O0 -Wall -g -fbacktrace -fbounds-check

#-- Intel Fortran --
F90 = ifort
OPT = -O3 -ip -ipo -unroll
#-no-inline-factor
#OPT = -O3
#OPT = -O0 -warn all -C -traceback
#OPT = -O0 -g

#OPT = -O3 -ip -ipo -arch SSE2 -tpp7 -xN # pablo.ethz.ch
#OPT = -O3 -ip -ipo -xB   # lapfabio
#OPT = -O3 -ip -ipo -unroll

#-- Lahey Fortran95 --
#F90 = lf95
#OPT = -O
#OPT = --chk aesux

#-- Digital --
#F90 = f95
#OPT = -fast

#== END USER SETTINGS =============================================

OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90)) 
#OBJS := $($(wildcard *.f90):.f90=.o) 

all: $(EXEC)

$(EXEC): $(OBJS)
	$(F90) $(OPT) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(EXEC) $(OBJS) *.M *.mod 

.SUFFIXES: 
.SUFFIXES: .f90 .o

%.o: %.f90
	$(F90) $(OPT) -c $<

include Makefile.depend
