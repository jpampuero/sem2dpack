%  SEM2DPACK/POST provides Matlab utilities for the manipulation
%  and visualization of SEM2DPACK simulations results.
%  
%  Reading simulation data:
%  
%   SEM2D_READ_SPECGRID	reads a spectral element grid
%   SEM2D_SNAPSHOT_READ	reads snapshot data
%   SEM2D_READ_SEIS 	reads seismogram data
%   SEM2D_READ_FAULT 	reads fault data
%  
%  Data manipulation:
%  
%   SEM2D_EXTRACT_POINT extracts field values at an arbitrary point
%   SEM2D_EXTRACT_LINE 	extracts field values along a vertical or horizontal line
%   ARIAS_INTENSITY 	computes Arias Intensity and Significant Duration
%   RESPONSE_SPECTRUM 	computes response spectra (peak dynamic response 
%   			of single-degree-of-freedom systems)
%  
%  Data visualization:
%  
%   SEM2D_PLOT_GRID 	plots a spectral element grid
%   SEM2D_SNAPSHOT_PLOT plots snapshot data
%   SEM2D_SNAPSHOT_GUI 	interactively plots snapshot data
%   SEM2D_SNAPSHOT_MOVIE makes an animation of snapshot data
%   PLOT_MODEL 		plots velocity and density model
%   PLOT_SEIS 		plots multiple seismograms
%   PLOT_FRONTS 	space-time plot of rupture front and process zone tail
%   SAMPLE_FAULT 	example of visualization of fault data
%   SAMPLE_SEIS 	example of visualization of seismogram data
%  
%  Miscellaneous tools:
%  
%   XCORRSHIFT 		cross-correlation time-delay measurement
%   SPECSHIFT 		signal time shift by non-integer lag via spectral domain
%   SPECFILTER 		zero-phase Butterworth filter via spectral domain
