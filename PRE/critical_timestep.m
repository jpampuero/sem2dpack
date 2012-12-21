% CRITICAL_TIMESTEP	Computes the critical timestep for the leapfrog time scheme
%
% SYNTAX	dtc = critical_timestep(csp,h,ngll)
%
% INPUT		csp	maximum wave speed (in m/s):
%			max S wave speed in SH mode (ndof=1) 
% 			max P wave speed in P-SV mode (ndof=2)
%		h	element size (in m)
%		ngll	number of GLL nodes per spectral element edge (2 to 20)
%			= polynomial order + 1
%
% OUTPUT	dtc	critical timestep
%
function dtc = critical_timestep(csp,h,ngll)

if ngll<2 | ngll>20, error('ngll must be from 2 to 20'); end

% tabulated values of critical frequency (non-dimensional, 1D)
% Omega_max(p) = Omega_max(ngll-1)
Omega_max = [2.0000000e+00 4.8989795e+00 8.6203822e+00 1.3540623e+01 1.9797952e+01 2.7378050e+01 ...
     3.6256848e+01 4.6421894e+01 5.7867306e+01 7.0590158e+01 8.4588883e+01 9.9862585e+01 ...
     1.1641072e+02 1.3423295e+02 1.5332903e+02 1.7369883e+02 1.9534221e+02 2.1825912e+02 2.4244948e+02];

% stability factor por leapfrog timescheme
C = 2;

% 2D critical time step, 
% assumes a square element
dtc = C*h/csp/sqrt(2)/Omega_max(ngll-1);
% in 3D replace sqrt(2) by sqrt(3)
