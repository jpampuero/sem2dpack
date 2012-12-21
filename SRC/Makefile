# SEM2DPACK version 2.2.12e -- A Spectral Element Method for 2D wave propagation and fracture dynamics,
#                             with emphasis on computational seismology and earthquake source dynamics.
# 
# Copyright (C) 2003-2007 Jean-Paul Ampuero
# All Rights Reserved
# 
# Jean-Paul Ampuero
# 
# ETH Zurich (Swiss Federal Institute of Technology)
# Institute of Geophysics
# Seismology and Geodynamics Group
# ETH Hönggerberg HPP O 13.1
# CH-8093 Zürich
# Switzerland
# 
# ampuero@erdw.ethz.ch
# +41 44 633 2197 (office)
# +41 44 633 1065 (fax)
# 
# http://www.sg.geophys.ethz.ch/geodynamics/ampuero/
# 
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
EXEC = $(HOME)/bin/sem2dsolve_12e

# 2. Uncomment/add/modify the lines corresponding to your compiler:

# GNU 
#F90 = g95
#OPT = 
#F90 = gfortran
#OPT = 

#-- Intel Fortran --
F90 = ifort
#OPT = -O3
#OPT = -O3 -ip -ipo -arch SSE2 -tpp7 -xN # pablo.ethz.ch
#OPT = -O3 -ip -ipo -xB   # lapfabio
OPT = -O3 -ip -ipo -unroll
#OPT = -warn all -C

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
