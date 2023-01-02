# Using sem2dpack to model site effects

## Highlights
* 2D sedimentary basins with irregular surface and/or sub-surface geometry  
* Soil nonlinearity of Iwan (XXX) 
* Cyclic mobility in liquefiable soils (Iai's liquefaction front model)
* Borehole boundary condition with incident waves
* Extra stations to write out nonlinearity features
* Post-processing jupyter notebooks including python scripts



## List of examples
### 1.1. Kathmandu Basin with obliquely incident plane waves

### 1.2. Kathmandu Basin with vertically incident plane waves


### 2. Hypothetical asymmetrical basin with liquefiable soil


### 3.1. In-plane (P-SV) Prenolin P1 benchmark model 
Single layer of nonlinear soil with borehole data from [PRENOLIN](xxx) benchmark  
Set the model parameters in the input file [`Par.inp`]
* In-plane model set by [`ndof=2`] in [`GENERAL`]
* [`IWAN`](XXX/Par.inp) material type
* Backbone curve set using spring number [`Nspr`] and reference strain [`gref`]
* Borehole condition set by [`borehole=T`] in [`BC_DIRNEU`]
* Input motion through borehole set in [`STF_TAB`]
* Extra station is set by [`extra=T`] in [`REC_LINE`] and detailed in [`REC_LINEX`]



### 3.2. Out-plane (SH) Prenolin P1 benchmark model 





## Tutorial steps
### 1. Code 
* Install and compile [`sem2dpack`](https://github.com/jpampuero/sem2dpack.git) with seismogenic width option


### 2. Example
* Set the model parameters in the input file [`Par.inp`](example_for_damage_2.5Dmodel/Par.inp)
* If wanted, modify station coordinates in [`stations`](example_for_damage_2.5Dmodel/stations)
* If wanted, modify model dimensions in [`layers`](example_for_damage_2.5Dmodel/layers)

### 3. Results
* To visulaise the simulation outputs, you can use the PY library in [`py-example`](py-example/)
* If necessary, modify the patth to simulation, in [`plot_fault_data`](py-example/plot_fault_data.py)
* Run the `plot_fault_data` to make plots of your choice. 



## Developer notes
* Each example was last tested in Caltech HPC cluster in December, 2022, 
with compiler option: xxx
* Only 'leap_frog' time scheme is tested



Oral, Weng, & Ampuero, 2019, Does a damaged fault zone mitigate the near-field
landslide risk during supershear earthquakes?—Application to the 2018 magnitude 7.5
Palu earthquake, Geophysical Research Letters, [DOI](https://doi.org/10.1029/2019GL085649), [PDF (preprint)](https://eartharxiv.org/repository/view/638/)

