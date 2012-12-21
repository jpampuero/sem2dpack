% PLOT_SEIS 	Plots multiple seismogram, verticallly offset
%
% SYNTAX	plot_seis(x,t,v,amp)
%		plot_seis(x,dt,v,amp)
%
% INPUT		x(nx)	linear position of receivers
%		t(nt)   time
%		dt	timestep
%		v(nt,nx) seismograms
%               amp     amplitude scaling factor = max amplitude / vertical offset
%                       (default=1)
%
function plot_seis(x,t,v,AFACTOR)

[NT,NX] = size(v);
if ~exist('AFACTOR','var') || isempty(AFACTOR), AFACTOR = 1; end

if length(x)~=NX, error('x and v have incompatible sizes'), end
if length(t)==1
  t = [0:NT-1]'*t;
elseif length(t)~=NT
  error('t and v have incompatible sizes')
end

ampli = 0.5*max(abs( v(:) ));
disp(sprintf('Amplitude scale (trace to trace) = %g',ampli))
if ampli>0
  offset = max(abs(diff(x)));
  ampli = AFACTOR*offset/ampli;
end
plot(t, v*ampli +repmat(x(:)',NT,1) );
title('SEM Seismograms')
xlabel('Time (s)')
ylabel('Distance (m)')
