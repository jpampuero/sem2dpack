# Using sem2dpack to model site effects

## Highlights
* 2D sedimentary basins with irregular surface and/or sub-surface geometry  
* Soil nonlinearity of Iwan (XXX) 
* Cyclic mobility in liquefiable soils (Iai's liquefaction front model)
* Borehole boundary condition with incident waves
* Extra stations to write out nonlinearity features
* Post-processing jupyter notebooks including python scripts in XXX

## List of examples and related publications
### 1. Prenolin P1 model 
Single layer of nonlinear soil with borehole data from [PRENOLIN](https://doi.org/10.1785/0120170210) benchmark  

E. Oral, 2016, Modélisation multi-dimensionnelle de la propagation des ondes sismiques dans des milieux linéaires et non-linéaires, Doctoral dissertation, Université Paris-Est. [PDF](https://theses.hal.science/tel-01562279/document)


## How-to set model properties in Par.inp file for each example
### 1.1. In-plane (P-SV) Prenolin P1 benchmark model 
* In-plane model by `ndof=2` in `GENERAL`
* `IWAN` material type
* Pressure-independent backbone curve prescribed by spring number `Nspr` and reference strain `gref`
* Borehole condition by `borehole=T` in `BC_DIRNEU`
* `hstf` and `vstf` for input in horizontal and vertical directions respectively
* Input motion(s) through borehole boundary set in `STF_TAB`  
* Extra stations set by `extra=T` in `REC_LINE` and detailed in `REC_LINEX`
### 1.2. Out-plane (SH) Prenolin P1 benchmark model 
Same as 1.1 example except for:
* `ndof=1` in `GENERAL`
* Only horizontal component for input motion

### 1.1. Kathmandu Basin with obliquely incident plane waves
### 1.2. Kathmandu Basin with vertically incident plane waves
### 2. Hypothetical asymmetrical basin with liquefiable soil



## Developer notes
* Each example was last tested in Caltech HPC cluster in December, 2022, 
with compiler option: xxx
* Only 'leap_frog' time scheme is tested

