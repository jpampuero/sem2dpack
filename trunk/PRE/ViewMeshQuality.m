% VIEWMESHQUALITY plots stability and resolution properties of a mesh
% 
% SYNTAX	[R,S]=ViewMeshQuality(0)
% 		[RI,SI]=ViewMeshQuality(1)
%
% INPUTS	plot_index	if =1 the stability and resolution measures 
%				are converted to logarithmic indices
%
% OUTPUTS	R(:)	a measure of the resolution of each element: 
%			an estimate of the local number of elements per wavelength (at 1 Hz) 
%			defined as
% 			R = minimum S wave speed / length of the largest element edge
%		RI(:)	a logarithmic "resolution index" defined for each element as
%  			RI = log10(R/median(R))
% 			where the median is taken over the whole mesh.
%		S(:)	a measure of the stability of each element:
% 			S = minimum (GLL node spacing / P wave velocity)
% 			The local critical timestep is proportional to S.
% 		SI(:)	a logarithmic "stability index" defined for each element as
%   			SI = log10(S/median(S))
% 			where the median is taken over the whole mesh.
%
% 		A figure is produced: 
%			The spatial distribution (left) and histogram (right)
%			of the stability measure S or it index SI (top) 
%			and the resolution measure R or its index RI (bottom)
%
function [Reso,Stab]=ViewMeshQuality(PLOT_INDEX);

if ~exist('PLOT_INDEX','var'), PLOT_INDEX = 1; end

% Load the finite element data
Knods = load('ElmtNodes_sem2d.tab');
Coorg = load('MeshNodesCoord_sem2d.tab');

% Load the resolution and stability check data 
Reso  = load('Resolution_sem2d.tab');
Stab  = load('Stability_sem2d.tab');
Stab = 1./Stab;

% convert to logarithmic index
if PLOT_INDEX
  Reso = log10( Reso/median(Reso) );
  Stab = log10( Stab/median(Stab) );
end

% Plot
map=colormap('jet'); map =map(end:-1:1,:); colormap(map);

subplot(2,3,[1 2])
  patch('Vertices',Coorg,'Faces',Knods(:,1:4), ...
        'FaceVertexCData',Stab,'FaceColor','flat', ...
        'EdgeColor','none');
axis equal
if PLOT_INDEX
  scale=max(abs(Stab));
  caxis([-1 1]*scale)
end
colorbar('SO')
title('Stability')

subplot(2,3,3)
hist(Stab)

subplot(2,3,[4 5])
  patch('Vertices',Coorg,'Faces',Knods(:,1:4), ...
        'FaceVertexCData',Reso,'FaceColor','flat', ...
        'EdgeColor','none');
axis equal
if PLOT_INDEX
  scale=max(abs(Reso));
  caxis([-1 1]*scale)
end
colorbar('SO')
title('Resolution')

subplot(2,3,6)
hist(Reso)
