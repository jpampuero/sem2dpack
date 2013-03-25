% SEM2D_VELOCITY_PLOT plots velocity at a given point from SEM2DPACK results
%
% SYNTAX     sem2d_velocity_plot(point)
%
% INPUT      point   location on fault for desired velocity plot
%                    default is 9.0km
%
%
function sem2d_velocity_plot(point)
clear d;
d = sem2d_read_fault();

%---Assumes point is at 9km------------
if ~exist('point','var')
   point = 9.0e3;
end
%--------------------------------------

%---Find spacing of grid along fault---
delta = zeros(length(d.x)-1,1);
for i=2:1:length(d.x)
    delta(i-1) = d.x(i)-d.x(i-1);
end
%--------------------------------------

%---Find index of input point----------
del=max(delta);
isel1=find(d.x > (point - del));
isel2=find(d.x < (point + del));
isel=size(intersect(isel1,isel2));
if isel(1,1)==3
    k=median(intersect(isel1,isel2));
elseif isel(1,1)==2
    k=min(intersect(isel1,isel2));
elseif isel(1,1)==1
    k=intersect(isel1,isel2);
else
    sprintf('Cannot locate %f along fault',point)
    return
end
%--------------------------------------

%---Plot velocity at point-------------
figure(1)
plot((0:d.nt-1)*d.dt,d.v(k,:))
title(sprintf('Velocity at x = %f km',(d.x(k)/1000)))
%--------------------------------------

end      


