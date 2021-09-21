% ************************************************************************
% This function generate the x and z coordinate of a rought fault

% Marion Thomas, last modified June 2018

%CALLS: 

%==========================================================================
%% ************************************************************************
function [x,z,n,s] = rough_flt_gen(rms,seed,BOX,h)

%% INPUTS VARIABLES

%Self-similar rought fault
% rms:      % root mean square (rms) height fluctuations of order "rms" times the profile length.
% nrms:     % name to assess the roughness (for the mesh grid name)
% seed: 	% Give a number or leave it empty for random generation of fault roughness
% BOX:      % domain limits in m [xmin xmax zmin zmax]
% h:       	% element size along fault in m (approximate)

%Local variables
maxiter=1e3;     	% Maximum iteration to find the right alpha
erriter=5;         	% Maximum difference (in %) bewteen the given and the generated alpha


%% FAULT GENERATION

%Number of elements
NELX = ceil(abs(BOX(2)-BOX(1))/h);      % number of elements in the x direction

%Define the self-similar rough fault
for n=1:maxiter
    [x,z,n,s,res] = oneDroughfault(h,NELX,rms,seed);
    if res.per < erriter
        disp('Fault Roughness')
        disp(' ')
        disp(['input alpha=',num2str(rms)]);  disp(['calculated alpha=',num2str(res.alpha)])
        disp(['ratio=',num2str(res.ratio)]);  disp([ 'difference(%)=', num2str(res.per)])
        break
    else
        [x,z,n,s,per] = oneDroughfault(h,NELX,rms,seed);
    end
end

