% ARIAS_INTENSITY computes Arias Intensity and Significant Duration
%
% SYNTAX	[Ia,SD] = arias_intensity(a,dt)
% 		[Ia,SD] = arias_intensity(a,dt,fig)
%
% INPUT		a(:)	accelerogram signal
%		dt	sampling timestep
%		fig	figure number to plot the cumulative Arias intensity as a function of time
%
% OUTPUT	Ia	Arias intensity = pi/(2g)*integral{a^2}
%		SD	Significant duration = time interval for the cumulative Ia(t)
%			to reach between 5% and 95% of its final value
%		Iat(:)	cumulative Arias intensity Ia(t)
%
function [Ia,SD,Iat] = arias_intensity(a,dt,fig)

g = 9.8; %acceleration of gravity
Iat = pi/(2*g)*dt*cumsum(a.^2);
Ia = Iat(end); % Arias intensity

% significant duration
% method 1, resolution = 1 sample
%[xmin,imin] = min( abs(Iat-0.05*Ia) );
%[xmax,imax] = min( abs(Iat-0.95*Ia) );
%SD = (imax-imin)*dt;
% method 2, subsample resolution
t = [0:length(Iat)-1]*dt;
  % remove constant portions, usually to handle leading zeros
  [Iat_bis,k]=unique(Iat,'last'); 
  t_bis = t(k);
tmin = interp1(Iat_bis,t_bis,0.05*Ia);
tmax = interp1(Iat_bis,t_bis,0.95*Ia);
SD = tmax-tmin;

% Husid plot Ia(t)
if exist('fig','var'),
  figure(fig)
  t = [0:length(Iat)-1]*dt;
 % left axis Arias intensity
  plot(t,Iat)
  xlabel('Time (s)')
  ylabel('Arias intensity (m/s)')
  c = axis;
  hold on
  plot([tmin tmin],[0 0.05*Ia],'k--',[c(1) tmin],[1 1]*0.05*Ia,'k--', ...
       [tmax tmax],[0 0.95*Ia],'k--',[c(1) tmax],[1 1]*0.95*Ia,'k--')
  hold off
 % right axis Ia(t)/Ia (%)
end
