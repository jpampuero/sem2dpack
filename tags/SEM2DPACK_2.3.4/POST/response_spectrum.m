% RESPONSE_SPECTRUM computes the peak dynamic response of single-degree-of-freedom systems
% (mass-spring-dashpot) using the Newmark Method for linear systems (Chopra's book page 177)
%
% SYNTAX	[maxA,maxV,maxD] = newm(a,dt,T,dr)
% 
% INPUT		a(:)	input accelerogram (in m/s^2)
%		dt	sampling timestep (in s)
%		T	fundamental period (in s)
%		dr	damping ratio = damping coefficient / critical damping coefficient
%
% OUTPUT	maxA	peak acceleration response 
% 		maxV	peak velocity response 
% 		maxD	peak displacement response 
%
% NOTE		The input "T" can be a vector. In that case the outputs are also vectors
%		(the response spectra).
%
% AUTHOR	Martin Mai (ETH Zurich)

function [maxA,maxV,maxD] = newm(ga,dt,T,dr);

% ------------------------------------------------------------------------------------------
% Newmark Method
% Variables initialization
  a = zeros(1,length(ga));
  v = zeros(1,length(ga));
  d = zeros(1,length(ga));
  
% Parameters defined for the Newmark method
% Linear acceleration method
  gamma = 1/2;
  beta = 1/4;

% Properties of the SDOF  
  w = 2*pi/T; % Circular frequency (rad/sec)
  m = 1; % mass
  k = m*w^2; % stiffness
  dc = 2 * dr * m * w; % damping coefficient (kN*s/m) (1 N = 1 kg-m/sec2)

% System starts from rest
% Initial calculations
  a(1) = (-m*ga(1) - dc*v(1) - k*d(1)) / m;
  
    kk = k + gamma * dc / beta / dt + m / beta / (dt^2); % variable k^ (kN/m)
  
     var_a = m / beta / dt + gamma * dc / beta; % variable a
     var_b = m / (2*beta) + dt * dc * (gamma/(2*beta) - 1); % variable b
  
% Calculation for each time step, i
  for j = 1 : (length(ga) - 1)
      dp = -m * (ga(j+1)-ga(j)) + var_a * v(j) + var_b * a(j); % kN
       du = dp / kk; % m
        dv = gamma * du / beta / dt - gamma * v(j) / beta + dt * a(j) * (1 - gamma/(2*beta)); % m/sec
         da = du / beta / (dt^2) - v(j) / beta / dt - a(j) / (2*beta); % m/sec2

      d(j+1) = d(j) + du;
      v(j+1) = v(j) + dv;
      a(j+1) = a(j) + da;
  end 

  maxA = max(abs(a+ga'));
  maxV = max(abs(v));
  maxD = max(abs(d));
  
return
