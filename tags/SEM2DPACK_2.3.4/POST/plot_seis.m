% PLOT_SEIS plots multiple seismogram traces, with a vertical offset
%
% SYNTAX	plot_seis(x,t,v)
%		plot_seis(x,dt,v)
%
% INPUT		x(nx)	linear position of receivers
%		t(nt)   time
%		dt	timestep
%		v(nt,nx) seismograms
%
function plot_seis(x,t,v)

[NT,NX] = size(v);

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
  ampli = offset/ampli;
end
plot(t, v*ampli +repmat(x(:)',NT,1) );
title('SEM Seismograms')
xlabel('Time (s)')
ylabel('Distance (m)')
