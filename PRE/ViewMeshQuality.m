% VIEWMESHQUALITY plots stability and resolution properties of a mesh

% Load the finite element data
Knods = load('ElmtNodes_sem2d.tab');
Coorg = load('MeshNodesCoord_sem2d.tab');
[NElem,ngnod] = size(Knods);
Ex=Coorg(Knods(:,1:4),1);
Ex=reshape(Ex,NElem,4)'; % +1.162e6; 
Ey=Coorg(Knods(:,1:4),2);
Ey=reshape(Ey,NElem,4)';

% Load the check data 
Reso  = load('Resolution_sem2d.tab');
Stab  = load('Stability_sem2d.tab');

% Modify to show relative values:
Reso = log10( Reso/median(Reso) );
Stab = -log10( Stab/median(Stab) );

% Plots
map=colormap('jet'); map =map(end:-1:1,:); colormap(map);

subplot(2,3,[1 2])
fill(Ex,Ey,Stab', 'LineStyle','none')
axis equal
scale=max(abs(Stab));
caxis([-1 1]*scale)
colorbar('SO')
title('Stability')

subplot(2,3,3)
hist(Stab)

subplot(2,3,[4 5])
fill(Ex,Ey,Reso', 'LineStyle','none')
axis equal
scale=max(abs(Reso));
caxis([-1 1]*scale)
colorbar('SO')
title('Resolution')

subplot(2,3,6)
hist(Reso)
