.SUFFIXES:
.SUFFIXES:  .f90 .f90.exe

BINPATH = $(HOME)/bin/
vpath %.exe $(BINPATH)

EXEC := $(patsubst %.f90,%.exe,$(wildcard *.f90)) 

F90 = ifort
OPT = -O3

all: $(EXEC)

%.exe: %.f90
	$(F90) $(OPT) $< -o $(BINPATH)$@ 

clean:
	cd $(BINPATH); rm -f $(EXEC)
