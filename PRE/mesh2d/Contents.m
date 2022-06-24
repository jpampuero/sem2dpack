% SEM2DPACK/PRE/MESH2D provides Matlab utilities for the generation, manipulation
% and visualization of structured 2D quadrilateral meshes, and unstructured
% compositions of structured meshes.
%
% The mesh database is stored in a structure described in MESH2D_TFI's help.
%
% Mesh generation:
%
%  MESH2D_TFI     	generates a structured mesh by transfinite interpolation
%  MESH2D_QUAD    	generates a structured mesh for a quadrilateral domain
%  MESH2D_CIRC_HOLE	generates a mesh for a square domain with a circular hole
%  MESH2D_WEDGE		generates a mesh for a triangular wedge domain
%  MESH2D_EX0		mesh for a vertical fault
%  MESH2D_EX1		mesh for a shallow layer over half-space with dipping fault
%  MESH2D_EX2   two horizontal faults cutting through a rectangular domain
%
% Mesh manipulation:
%
%  MESH2D_ROTATE   	rotates the node coordinates
%  MESH2D_TRANSLATE	translates the node coordinates
%  MESH2D_MERGE    	merges several meshes into a single mesh
%  MESH2D_WRITE  	writes a 2D mesh database file (*.mesh2d)
%  MESH2D_READ 		reads a 2D mesh database from a *.mesh2d file
%
% Reading mesh files from other mesh generation software:
%
%  READ_DCM 		reads a 2D mesh in the DCM format of EZ4U (http://www-lacan.upc.es/ez4u.htm)
%  READ_INP		reads a 2D mesh in the INP format of ABAQUS exported by CUBIT
%
% Mesh visualization:
%
%  MESH2D_PLOT        	plots a 2D mesh
%
% Miscellaneous tools:
%
%  SAMPLE_SEGMENTS	generates points that regularly sample multiple segments of a line
%
