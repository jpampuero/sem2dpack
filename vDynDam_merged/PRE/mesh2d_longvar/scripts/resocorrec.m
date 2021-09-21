% ************************************************************************
% This function generate the x and z coordinate of a rought fault

% Marion Thomas, last modified June 2018

%CALLS: 

%==========================================================================
%% ************************************************************************
function [x] = resocorrec(x,h,Q)

%% INPUTS VARIABLES

%x      % variables to be corrected
% h     % element size in m (resolution)
% Q     % quad element type, identified by its number of nodes: 4 or 9.


%% FAULT GENERATION

resadj=Q*h;
x=ceil(x./resadj).*resadj;

