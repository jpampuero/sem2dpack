% SET_MESH_LAYERS sets the parameters of a layered mesh for 2D SH waves
% according to a minimum number of elements per wavelength
%
% INPUTS	model	layered media structure with the following components
%			for each layer, listed from top to bottom:
%			rho	density
%			cs	S wave speed
%			h	layer thickness
% 		fmax 	maximum frequency to be resolved
% 		Npw	[1] minimum number of elements per wavelength
%		fig	[1] figure number (0 = no plot)
%
% OUTPUTS	p	mesh parameters structure, containing:
%			xlim(2)	limits of the domain along X
%			zlim(2)	limits of the domain along Z
%			nx	number of elements along X
%			nz(:)	number of elements along Z for each layer
%			ztop(:) Z position of the top of each layer
%			tag(:)	material tag for each layer
%
function p = set_mesh_layers(model,fmax,Npw,fig)

VERSION =1; % 0= MESH_LAYER blocks, 1= separate file

if ~exist('Npw','var'), Npw=1; end
if ~exist('fig','var'), fig=1; end

nl = length(model.cs);
LX = model.xlim(2)-model.xlim(1);
lmin = model.cs/fmax;

p.ztop = [0; -cumsum(model.h(1:end-1))];
p.xlim = model.xlim;
p.zlim = [0 -sum(model.h)];
p.tag = [1:nl]';
p.nx = ceil( max(Npw*LX./lmin) ); 
p.nz = ceil(Npw*model.h./lmin);

%------ export input files ------

disp('&MESH_DEF method = "LAYERED" /')
if VERSION==0
  disp(sprintf('\n&MESH_LAYERED xlim=%f,%f, zmin=%f, nx=%u, nlayer=%u /\n', ...
       p.xlim, p.zlim(2), p.nx, nl))
  disp(sprintf('&MESH_LAYER nz=%u, ztop=%f, tag=%u /\n', [p.nz,p.ztop,p.tag]'))     

else
  disp(sprintf('\n&MESH_LAYERED xlim=%f,%f, zmin=%f, nx=%u, file="%s" /\n', ...
       p.xlim, p.zlim(2), p.nx, 'layers.tab'))
  fid=fopen('layers.tab','w');
  fprintf(fid,'%f %u %u\n', [p.ztop,p.nz,p.tag]');
  fclose(fid);
end

disp(sprintf('&MATERIAL tag=%u, kind=''ELAST'' /\n &MAT_ELASTIC rho=%f, cp=%f, cs=%f /\n', ...
     [p.tag(:),model.rho(:),model.cp(:),model.cs(:)]'))     

%------ plots --------
if ~fig, return, end

cs2=[model.cs(:) model.cs(:)]'; 
cs2=cs2(:);

h = model.h;
h(end) = cs2(end)/fmax;
h2=[h' ;h'];
h2=h2(:);

dep=cumsum(model.h);
dep2=[[0;dep(1:end-1)]';dep'];  
dep2=dep2(:);

hol = h2./(cs2/fmax);
nel2 = ceil(Npw*hol);

dtx = LX/p.nx ./cs2;
dtz = h2./nel2./cs2;
dtz = dtz;
dt = max([dtx(:);dtz(:)]);
dtx = dtx/dt;
dtz = dtz/dt;

figure(fig)

subplot(131)
plot(cs2,-dep2)
ylabel('z (m)')
xlabel('c_s (m/s)')
title('S wave velocity model')
axis([0 1500 -inf 0])

subplot(132)
plot( nel2,-dep2, '--',hol,-dep2 , nel2./hol, -dep2, ':')
xlabel('H / \lambda_{min} = H f_{max}/c_s')
legend('Elements per layer','layer thickness / \lambda_{min}','Elements per \lambda_{min}',4)
title(sprintf('%s\n%s',...
     'For mesh generation with f_{max} = 20 Hz',...
     '(resolution criterion)'))
axis([0 7 -inf 0])


subplot(133)
plot( dtx, -dep2, dtz,-dep2)
xlabel('\Delta t_c / max( \Delta t_c )')
title(sprintf('%s\n%s',...
     'For time step selection',...
     '(stability criterion)'))
legend('X','Z',4)
axis([0 1.2 -inf 0])

