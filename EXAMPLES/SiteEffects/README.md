# Using sem2dpack to model site effects

## New features
* Examples for 2D sedimentary basins with irregular surface and/or sub-surface geometry  
* Soil nonlinearity of [Iwan (1967)](https://doi.org/10.1115/1.3607751)
* Viscoelasticity of [Liu and Archuleta (2006)](https://doi.org/10.1785/0120050173)
* Cyclic mobility in liquefiable soils by the liquefaction front model of [Iai et al. (1992)](https://doi.org/10.3208/sandf1972.32.2_1)
* Borehole boundary condition with incident waves
* Extra stations to write out soil nonlinearity features (stress-strain curves, number of activated plasticity surfaces)
* Post-processing with jupyter notebooks including python scripts in [sem2dpack/JUPYTER](https://github.com/elifo/sem2dpack/tree/master/JUPYTER/SiteEffects)


## List of examples and related publications
### 1. Prenolin P1 model 
Single layer of nonlinear soil with borehole data from [PRENOLIN](https://doi.org/10.1785/0120170210) benchmark  
  * E. Oral, 2016, Modélisation multi-dimensionnelle de la propagation des ondes sismiques dans des milieux linéaires et non-linéaires, Doctoral dissertation, Université Paris-Est. [PDF](https://theses.hal.science/tel-01562279/document)

### 2. Hypothetical asymmetrical basin model of Oral et al. (2019)
2D multi-layered model with an assymetrical shape, useful to approximate Earth structures such as Alpine valleys and sedimentary basins in New Zealand 
  * E. Oral, C. Gélis, and L.F. Bonilla, 2019, 2-D P-SV and SH spectral element modelling of seismic wave propagation in non-linear media with pore-pressure effects, Geophysical Journal International, 217(2), pp.1353-1365. [DOI](https://doi.org/10.1093/gji/ggz041)
  * E. Oral, 2016, Modélisation multi-dimensionnelle de la propagation des ondes sismiques dans des milieux linéaires et non-linéaires, Doctoral dissertation, Université Paris-Est. [PDF](https://theses.hal.science/tel-01562279/document)

### 3. Kathmandu basin model of Oral et al. (2022)
2D multi-layered model of Kathmandu City, Nepal, where site effects including soil nonlinearity were observed during the 2015 earthquake. 
  * E. Oral, P. Ayoubi, J.P. Ampuero, D. Asimaki, and L.F. Bonilla, 2022, Kathmandu Basin as a local modulator of seismic waves: 2-D modelling of non-linear site response under obliquely incident waves, Geophysical Journal International, 231(3), pp.1996-2008. [DOI](https://doi.org/10.1093/gji/ggac302)


## Model properties set in Par.inp file for each example

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

### 2. Hypothetical asymmetrical basin with liquefiable soil
* In-plane model by `ndof=2` in `GENERAL`
* Viscoelasticity set by `VEPMOD=T` in `MAT_IWAN` for nonlinear soils (inside basin)
* Pressure-dependent nonlinearity model with failure line angle `phi_f` provided
* Viscoelasticity parameters `Qp`, `Qs`, and `fr` (quality factors and reference frequency) given in `MAT_VEP`
* Water table level set by once `WT` in `MAT_IWAN`
* Liquefiable top layer (tag=1) set by `IAIMOD=T` in `MAT_IWAN`
* Iai model parameters given in `MAT_IAI`
* `VISLA` material type for no nonlinearity (outside basin) with model parameters given in `MAT_VISLA`
* Extra stations set by `extra=T` in `REC_LINE` and detailed in `REC_LINEX`
* Lateral boundaries are set to periodic in `BC_DEF`
* Vertical plane wave incidence through absorbing bottom boundary in `BC_ABSORB`
* Incident motion set in `SRC_WAVE` and `STF_TAB`

### 3.1. Kathmandu Basin with vertically incident plane waves
* In-plane model by `ndof=2` in `GENERAL`
* Nonlinear soil of top layer (tag=1) adjusted to Nepal topography by `IS_NEPAL=T` in `MAT_IWAN`
* Pressure-dependent model with failure line angle `phi_f` and `cohesion` given
* `VISLA` material type for all other layers with model parameters given in `MAT_VISLA`
* Extra stations set by `extra=T` in `REC_LINE` and detailed in `REC_LINEX`
* Vertical plane wave incidence through absorbing layer in bottom in `BC_ABSORB`
* Incident motion set in `SRC_WAVE` and `STF_TAB`

### 3.2. Kathmandu Basin with obliquely incident plane waves
Same as 3.1 example except for:
* Lateral boundaries set to absorbing (tags=2,4) in `BC_ABSORB`
* Oblique plane wave incidence `Angle = 30` in `SRC_WAVE`


## Notes for developers/advanced users
* Each example was tested in Caltech HPC cluster in 12/2022 with compilation: ifort -O3
* Only `leap_frog` time scheme is used and tested for highlighted new features


fin
